##################################################################
# Sieve Estimator For VE waning
##################################################################

library(tidyverse)
library(mgcv)
library(splines)
library(coneproj)

psi_bs_early <- function(Tvec, Vvec, Vvec_early, df, left_knot = 0){
  
  K = df - 3
  
  cuts <- seq(1, K) / (K + 1)
  
  knots <- quantile(pmax(0, Tvec-Vvec)[which(Tvec > Vvec)], cuts)
  
  boundary.knots = c(left_knot, max(pmax(0, Tvec-Vvec)))
  
  basis <- bs(pmax(0, Tvec-Vvec), df=df, degree = 3, knots=knots, Boundary.knots = boundary.knots, intercept=TRUE)
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      as.numeric(Tvec - Vvec > 0) * as.matrix(basis)),
    basis = basis
  ))
  
}

sieve_multinomial_EM <- function(dat, J.name="J", V.name="V", V.early.name=NULL, T.name="T", 
                              psi_delta, verbose=TRUE, monotone=FALSE, Amat=NULL, maxiter=50, tol=1E-6, calc_mixture = TRUE, pk_names = NULL, ...){
  
  # J.name="J" ; T.name="T" ; V.name="V"; V.early.name="V_early"
  # 
  # dat=events %>% dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
  #   mutate(J = case_when(
  #     type=="Omicron" ~ 2,
  #     type=="Delta" ~ 1,
  #     type=="Non-SARS" ~ 0
  #   )) %>%
  #   rename(`T`=tstop, V=Xstart_full, V_early = Xstart_early)
  
  #######
  # Create basis for f()
  #######
  
  n <- nrow(dat)
  J <- as.numeric(dat[[J.name]])
  Tvec <- dat[[T.name]]
  Vvec <- dat[[V.name]]
  Vvec <- ifelse(is.na(Vvec), Inf, Vvec)
  # build Z matrix (nxd)
  if(!is.null(V.early.name)){
    Vvec_early <- dat[[V.early.name]]
    Vvec_early <- ifelse(is.na(Vvec_early), Inf, Vvec_early)
    base<-psi_delta(Tvec, Vvec, Vvec_early, ...)
  }else{
    base<-psi_delta(Tvec, Vvec, ...)
  }
  
  #base<-psi_bs_early(Tvec, Vvec, Vvec_early, df=5)
  Z <- base$mat
  d <- ncol(Z)
  colnames(Z)<-paste0("f", 1:d)
  
  #######
  # Create basis for alpha
  #######
  
  # Build basis matrix for alpha
  K = floor(n^(1/3))
  alpha_basis <- splines::bs(dat[[T.name]], df = K+2, degree=3)
  alphas<-as.data.frame(predict(alpha_basis, dat[[T.name]]))
  #alphas <- apply(alphas, MARGIN=2, FUN=function(x){x - mean(x)})
  colnames(alphas)<-paste0("alpha", 1:ncol(alphas))
  
  #######
  # Create multinomial log likelihood function
  #######
  
  loglik <- function(params, J, alpha, Z, offset_1, offset_2) {
    
    #params <- fit$par
    #params<-init
    p <- ncol(Z)
    q <- ncol(alpha) #number of shared alpha columns
    beta1 <- params[1:p]
    beta2 <- params[(p+1):(2*p)]
    gamma <- params[((2*p)+1):length(params)]
    
    alpha <- as.matrix(alpha)
    
    #dim(alpha)
    #dim(gamma)
    
    # Linear predictors
    eta1 <- as.numeric(alpha %*% gamma + Z %*% beta1)
    eta2 <- as.numeric(alpha %*% gamma + Z %*% beta2)
    
    # Probabilities
    denom <- 1 + (offset_1 * exp(eta1)) + (offset_2 * exp(eta2))
    p0 <- 1/denom
    p1 <- (offset_1 * exp(eta1))/denom
    p2 <- (offset_2 * exp(eta2))/denom
    
    # Log-likelihood
    ll <- sum(log(ifelse(J==0, p0, ifelse(J==1, p1, p2))))
    return(-ll)  # negative for minimization
  }
  
  conv_flag = FALSE
  iter=0
  while(iter<=maxiter & conv_flag==FALSE){
  
  iter=iter+1
    
  print(iter)
  #Iterate between estimating variant mixture proportions and model output
    
  if(iter==1){
    
    if(calc_mixture==TRUE){
    to_regress_p <- bind_cols("T"=dat[[T.name]], 
                              "J"=dat[[J.name]]) %>% 
      bind_cols(., as.data.frame(Z)) %>%
      filter(J != 0) %>%
      mutate(J = J-1) %>%
      mutate(m_col_sum = across(starts_with("f")) %>% rowSums()) %>%
      filter(m_col_sum==0)
    
    fit_p<-gam(J ~ s(`T`), data=to_regress_p, family=binomial)
    
    pk <- predict(fit_p, newdata = bind_cols("T"=dat[[T.name]]), type="response")
    }else{
    pk <- dat[[pk_names]]
    }
    
    init <- rep(0, ncol(alphas) + 2*ncol(Z))  # gamma + beta1 + beta2
  
    #######
    # Create optimize for parameter vector
    #######
    
    fit <- optim(init, 
               loglik, 
               J=J, 
               alpha=as.matrix(alphas), 
               Z=as.matrix(Z %>% as.data.frame() %>% dplyr::select(starts_with("f"))), 
               offset_1 = (1-pk), 
               offset_2 = (pk), 
               method="BFGS", 
               hessian=TRUE)
    
    beta <- fit$par[1:(2*ncol(Z))]
    beta_1 <- fit$par[1:(ncol(Z))]
    beta_2 <- fit$par[(ncol(Z)+1):(2*ncol(Z))]
    next
  }
  
  if(iter!=1){
    
  if(calc_mixture==TRUE){
   a_k = pk * exp(as.numeric(Z %*% beta_2)) / (pk * exp(as.numeric(Z %*% beta_2)) + (1-pk) * exp(as.numeric(Z %*% beta_1)))
   
   fit <- optim(init, 
                loglik, 
                J=J, 
                alpha=as.matrix(alphas), 
                Z=as.matrix(Z %>% as.data.frame() %>% dplyr::select(starts_with("f"))),, 
                offset_1 = (1-a_k), 
                offset_2 = (a_k), 
                method="BFGS", 
                hessian=TRUE)
   
   beta_new <- fit$par[1:(2*ncol(Z))]
   beta_1_new <- fit$par[1:(ncol(Z))]
   beta_2_new <- fit$par[(ncol(Z)+1):(2*ncol(Z))]
   
   fit_p<-gam(a_k ~ s(Tvec), data=NULL, family=binomial)
   
   pk_new = predict(fit_p, newdata=data.frame(Tvec), type="response")
  }else{
    
    fit <- optim(init, 
                 loglik, 
                 J=J, 
                 alpha=as.matrix(alphas), 
                 Z=as.matrix(Z %>% as.data.frame() %>% dplyr::select(starts_with("f"))),, 
                 offset_1 = (1-pk), 
                 offset_2 = (pk), 
                 method="BFGS", 
                 hessian=TRUE)
    
    beta_new <- fit$par[1:(2*ncol(Z))]
    beta_1_new <- fit$par[1:(ncol(Z))]
    beta_2_new <- fit$par[(ncol(Z)+1):(2*ncol(Z))]
  
    pk_new <- dat[[pk_names]]
  }
  }
    
  print(median(sqrt((beta-beta_new)^2)))
  print(mean((pk-pk_new)^2))
  conv_flag = (median(sqrt((beta-beta_new)^2)) <= tol) & (mean((pk-pk_new)^2) <= tol)
    
  beta = beta_new
  beta_1 = beta_1_new
  beta_2 = beta_2_new
  pk = pk_new
  }

  vcov <- solve(fit$hessian)
  
  se_beta <- sqrt((diag(vcov)[1:(2*ncol(Z))]))
  
  if(monotone==FALSE){
    Amat=NULL
    beta_mono=vector(length=0)
  }
  
  if(monotone==TRUE){
    if (verbose) cat(sprintf("Smoothing beta to be monotone"))
    
    beta_mono_1 <- as.numeric(coneproj::coneA(beta[1:(ncol(Z))], amat = Amat, w = (1/(diag(vcov)[1:ncol(Z)])))$thetahat)
    
    beta_mono_2 <- as.numeric(coneproj::coneA(beta[(ncol(Z)+1):(2*ncol(Z))], amat = Amat, w = (1/(diag(vcov)[(ncol(Z)+1):(2*ncol(Z))])))$thetahat)
    
  }
  
  if(monotone==TRUE){
    
    beta_mono_1 = beta_mono_1
    beta_mono_2 = beta_mono_2
    
  }else{
    
    beta_mono_1 = beta_1
    beta_mono_2 = beta_2
    
  }
  
  
  return(
    list(
      beta_1 = beta_mono_1,
      beta_2 = beta_mono_2,
      beta_1_unconstr = beta_1,
      beta_2_unconstr = beta_2,
      se = se_beta,
      cov = vcov,
      Amat = Amat,
      basis = base$basis
    )
  )
}

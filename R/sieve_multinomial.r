#' partially_linear_multinomial_sieve
#'
#' Function used to estimate time-varying VE using the partially linear *multinomial* logistic regression model adjusting for vaccine-irrelevant infections via the sieve estimation method.
#' 
#' The code has an option to profile out nuisance parameters (the variant mixture proportions) by using external disease surveillance data.
#' 
#' @param dat data.frame, containing columns J.name, V.name, T.name, and V.early.name (optional).
#' @param J.name character, column name associated with outcome type (J=0 vaccine-irrelevant, J=1 vaccine-preventable)
#' @param V.name character, column name associated with vaccination date or peak immunization date (if V.early.name is specified)
#' @param V.early.name character, OPTIONAL column name associated with "early" effectiveness (often the vaccination date itself)
#' @param T.name character, column name associated with infection date
#' @param S.name character, column name associated with the strain name (S=NA for irrelevant, S=1 for strain 1, S=2 for strain 2, etc.)
#' @param psi_delta function, used to build the basis for time-varying VE
#' @param m integer, number of strains of vaccine-preventable infection in the dataset
#' @param p_s_fun function, OPTIONAL function that takes as input a vector of times and outputs mixture probabilities for the different strains
#' @param calc_mixture logical, whether or not to calculate variant mixture proportions from the data itself or rely on external source
#' @param tol numeric, tolerance limit to stop iterations
#' @param verbose logical, whether to print out intermediate info
#' @param ... extra arguments to use for basis creation functions
#' @return list containing parameter estimates associated with VE, SEs, vcov matrix, the varying intercept alpha(t), and the basis
#' @export

partially_linear_multinomial_sieve <- function(dat, 
                                               J.name="J", 
                                               V.name="V", 
                                               V.early.name=NULL, 
                                               T.name="T", 
                                               S.name="S",
                                               m,
                                               psi_delta, 
                                               p_s_fun = NULL,
                                               verbose=TRUE, 
                                               tol=1E-6, 
                                               calc_mixture = FALSE, ...){
  
  n <- nrow(dat) # number of observations
  J <- dat[[J.name]] # Vaccine preventable/irrelevant event labels
  Tvec <- dat[[T.name]] # Event times
  Vvec <- dat[[V.name]] # Vaccination times
  Svec <- dat[[S.name]] # Strain labels (0 if off-target)
  
  m=m
  
  bases<-vector(mode="list", length=m)
  
  for(i in 1:m){
    if(!is.null(V.early.name)){
      Vvec_early <- dat[[V.early.name]]
      bases[[i]]<-psi_delta(Tvec, Vvec, Vvec_early, ...)
    }else{
      bases[[i]]<-psi_delta(Tvec, Vvec, ...)
    }
  }
  
  Z_stack <- do.call(cbind, lapply(1:length(bases),
                                   FUN=function(x){
                                     bases[[x]]$mat
                                   }))
  
  p_tot <- ncol(Z_stack)
  d = p_tot/m
  
  #######
  # Create basis for alpha
  #######
  
  # Build basis matrix for alpha
  K = floor(n^(1/3))
  alpha_basis <- splines::bs(Tvec, df = K+2, degree=3)
  alphas<-as.data.frame(predict(alpha_basis, Tvec))
  colnames(alphas)<-paste0("alpha", 1:ncol(alphas))
  
  #Initialize multinomial log likelihood function
  
  loglik <- function(params, J, S, m, bases, alpha, offsets) {
    
    alpha <- as.matrix(alpha)
    p <- ncol(bases$mat[[1]])
    q <- ncol(alpha) #number of shared alpha columns

    gamma <- params[(p*m + 1):length(params)]
    exp_eta <- matrix(0, nrow=n, ncol=m+1)
    
    for(i in 1:m){
      
    Z<-bases$mat[[i]]
    beta_s <- params[(p*(i-1)+1):(p*(i-1)+p)]
    
    eta <-as.numeric(alpha %*% gamma + Z %*% beta_s)
    
    exp_eta[,i] = offsets[,i] * exp(eta)
    
    }
    
    exp_eta[,m+1] <- 1
    
    denoms <- rowSums(exp_eta)
    
    #matrix with P(J=1, S=1), P(J=1, S=2), ..., P(J=0)
    probs <- exp_eta/denoms
    
    ll_terms<-rep(0, nrow(probs))
    for(i in 1:length(J)){
      
      ll_terms[i] <- dplyr::case_when(
        J[i] == 0 ~ log(probs[i,(m+1)]),
        J[i] == 1 ~ log(probs[i, S[i]])
      )
      
    }

    # Log-likelihood
    ll <- sum(ll_terms)
    return(-ll)  # negative for minimization
  }
  
  # Whether or not to calculate mixture proportions based on study data OR to rely on external data
  
  if(calc_mixture==TRUE){
    
    if(is.null(V.early.name)){
    to_regress_p <- dplyr::bind_cols("T"=dat[[T.name]], 
                              "J"=dat[[J.name]],
                              "S"=dat[[S.name]],
                              "V"=dat[[V.name]]) %>%
      dplyr::filter(`T` <= V)
    
    fit_p<-mgcv::gam(S ~ s(`T`), data=to_regress_p, family=mgcv::multinom(K=m))
    }else{
    to_regress_p <- dplyr::bind_cols("T"=dat[[T.name]], 
                                "J"=dat[[J.name]],
                                "S"=dat[[S.name]],
                                "V"=dat[[V.early.name]]) %>%
      dplyr::filter(`T` <= V)
    
    fit_p<-mgcv::gam(S ~ s(`T`), data=to_regress_p, family=mgcv::multinom(K=m))
    }
    
    pk <- predict(fit_p, newdata = dplyr::bind_cols("T"=dat[[T.name]]), type="response")
  }else{
    pk <- p_s_fun(Tvec)
  }
  
  init <- rep(0, ncol(alphas) + ncol(Z_stack))  # gamma + beta1 + beta2
  
  # Find parameters that maximize log likelihood
  
  fit <- stats::optim(init, 
                      loglik, 
                      J=J,
                      S=S,
                      m=m,
                      bases = bases,
                      alpha=as.matrix(alphas), 
                      offsets = pk, 
                      control = list(abstol=tol),
                      method="BFGS", 
                      hessian=TRUE)
  
  p<-ncol(bases$mat[[1]])
  betas <- vector(mode="list", length=m)
  vcov <- solve(fit$hessian)
  se_beta <- sqrt((diag(vcov)[1:(m*p)]))
  for(i in 1:m){
    betas[[i]] <- params[(p*(i-1)+1):(p*(i-1)+p)]
  }
  
  return(
    list(
      betas=betas,
      se = se_beta,
      cov = vcov,
      bases = bases
    )
  )
}

##################################################################
# Sieve Estimator For VE waning
##################################################################

# library(tidyverse)
# library(mgcv)
# library(splines)
# library(coneproj)

# psi_d2 <- function(Tvec, Vvec) {
#   # Vmat: n x d
#   # default: psi(T-V) = V
#   return(cbind(as.numeric(Tvec - Vvec > 0), pmax(0, Tvec-Vvec)))
# }
# 
# psi_bs <- function(Tvec, Vvec, df=3) {
#   # Vmat: n x d
#   # default: psi(T-V) = V
#   return(cbind(as.numeric(Tvec - Vvec > 0), bs(pmax(0, Tvec-Vvec), df=df, degree = 3)))
# }
# 
# psi_d2_early <- function(Tvec, Vvec, Vvec_early){
#   
#   return(cbind(
#     as.numeric(between(Tvec, Vvec_early, Vvec)), #early vaccination
#     as.numeric(Tvec - Vvec > 0), #full vaccination
#     pmax(0, Tvec-Vvec))) #time-since-full vaccination
#   
# }

sieve_partially_linear_logistic <- function(dat, J.name="J", V.name="V", V.early.name=NULL, T.name="T", 
                                            psi_delta, verbose=TRUE, monotone=FALSE, Amat=NULL, ...){
  
  
  #dat<- dat_const_run

  # dat=delta_eligible
  #J.name="J" ; T.name="T" ; V.name="V"; V.early.name="V_early"
  
  n <- nrow(dat)
  J <- as.numeric(dat[[J.name]])
  Tvec <- dat[[T.name]]
  Vvec <- dat[[V.name]]
  # build Z matrix (nxd)
  if(!is.null(V.early.name)){
    Vvec_early <- dat[[V.early.name]]
    base<-psi_delta(Tvec, Vvec, Vvec_early, ...)
  }else{
    base<-psi_delta(Tvec, Vvec, ...)
  }
  
  Z <- base$mat
  
  d <- ncol(Z)
  colnames(Z)<-paste0("f", 1:d)
  
  # Build basis matrix for alpha
  K = floor(n^(1/3.5))
  alpha_basis <- splines::bs(dat[[T.name]], df = K+2, degree=3)
  
  alphas<-as.data.frame(predict(alpha_basis, dat[[T.name]]))
  #alphas <- apply(alphas, MARGIN=2, FUN=function(x){x - mean(x)})
  colnames(alphas)<-paste0("alpha", 1:ncol(alphas))
  
  to_regress <- bind_cols("J"=dat[[J.name]], Z, alphas)
  
  fit<-glm(as.formula(
  paste0("J ~ -1 + ", paste0("f", 1:ncol(Z), collapse="+"), "+", paste0("alpha", 1:ncol(alphas), collapse="+"), collapse="")),
  data=to_regress, family=binomial)
  
  beta <- coef(fit)[1:d]
  se_beta <- summary(fit)$coef[1:d,2]
  cov <- vcov(fit)[1:d, 1:d]
  alpha <- coef(fit)[(d+1):length(coef(fit))]
  
  # update: Smooth beta
  if(monotone==FALSE){
    Amat=NULL
    beta_mono=vector(length=0)
  }
  if(monotone==TRUE){
    if (verbose) cat(sprintf("Smoothing beta to be monotone"))
    
    
    if(is.null(Amat)){
      Amat <- matrix(0, nrow=length(beta), ncol=length(beta))
      
      # This works for B-spline *without* I(V ≤ T) intercept
      for(i in 1:nrow(Amat)){
        
        if(length(beta)==2){
          
          Amat[1,1] = -1
          Amat[2,2] = 1
          
        }else{
          
          if(i==1){
            Amat[i,i] = -1
          }else{
            Amat[i,(i-1)] = -1
            Amat[i,(i)] = 1
          }
        }
      }
    }else{
      Amat = Amat
    }
    
    beta_mono <- coneproj::coneA(y=beta, amat=Amat, w=1/se_beta^2)$theta[,1]
  }
  
  if(monotone==TRUE){
    
    beta_mono = beta_mono
    
  }else{
    
    beta_mono=beta
    
  }
  
  return(
    list(
      beta = beta_mono,
      beta_unconstr = beta,
      se = se_beta,
      cov = cov,
      Amat = Amat,
      basis = base$basis,
      alpha = alpha
    )
  )
  
}

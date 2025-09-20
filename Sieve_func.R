##################################################################
# Sieve Estimator For VE waning
##################################################################

library(tidyverse)
library(mgcv)
library(splines)
library(coneproj)

psi_d2 <- function(Tvec, Vvec) {
  # Vmat: n x d
  # default: psi(T-V) = V
  return(cbind(as.numeric(Tvec - Vvec > 0), pmax(0, Tvec-Vvec)))
}

psi_bs <- function(Tvec, Vvec, df=3) {
  # Vmat: n x d
  # default: psi(T-V) = V
  return(cbind(as.numeric(Tvec - Vvec > 0), bs(pmax(0, Tvec-Vvec), df=df, degree = 3)))
}

sieve_partially_linear_logistic <- function(dat, J.name="J", V.name="V", T.name="T",
                                            psi_delta, add_denom_ridge = 1e-8, verbose=TRUE, monotone=FALSE){
  
  
  
  # dat <- dat_const %>%
  #   filter(Delta_Y==1 | Delta_N==1) %>%
  #   transmute(
  #     `T` = pmin(Y, N),
  #     `A` = as.numeric(Z <= T),
  #     `J` = ifelse(`T`==Y, 1, 0),
  #     `V` = Z
  #   )
  #dat <- dat_wane_run
  
  n <- nrow(dat)
  
  # Build psi for VE waning
  
  psi <- psi_delta(dat[[T.name]], dat[[V.name]])
  d <- ncol(psi)
  colnames(psi) <- paste0("f", 1:d)
  Zmat <- as.data.frame(psi)
  
  # Build 
  K = floor(n^(1/4))
  alpha_basis <- splines::bs(dat[[T.name]], df = K+2, degree=3)
  
  alphas<-as.data.frame(predict(alpha_basis, dat[[T.name]]))
  alphas <- apply(alphas, MARGIN=2, FUN=function(x){x - mean(x)})
  colnames(alphas)<-paste0("alpha", 1:ncol(alphas))
  
  to_regress <- bind_cols("J"=dat[[J.name]], Zmat, alphas)
  
  fit<-glm(as.formula(
  paste0("J ~ ", paste0("f", 1:ncol(Zmat), collapse="+"), "+", paste0("alpha", 1:ncol(alphas), collapse="+"), collapse="")),
  data=to_regress, family=binomial)
  
  if(monotone==FALSE){
  
  beta <- coef(fit)[2:(d+1)]
  se_beta <- summary(fit)$coef[2:(d+1),2]
  cov <- vcov(fit)[2:(d+1), 2:(d+1)]
  
  list(beta = beta,
       se = se_beta,
       cov = var_beta)
  }
  
  if(monotone==TRUE){
  
  amat<-matrix(nrow=length(coef(fit))-1, ncol=length(coef(fit)))
  for(i in 1:nrow(amat)){
    
    if(i<=2){
      
      amat[i,]<-0
      
    }else if(i>2 & i <= d){
      
      amat[i,]<-c(rep(0, i-1), -1, 1, rep(0, ncol(amat)-(i+1)))
      
    }else{
      
      amat[i,]<-0
    }
    
  }
  
  coneproj::coneA(coef(fit), amat=amat)
  }
  
}
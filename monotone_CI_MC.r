##########
# Depreciated method for Monte Carlo sampling from asymptotic dist of unrestricted estimator

monotone_CI_MC <- function(beta, vcov, Amat, w, grid, M=5000, seed=47){
  
  set.seed(seed)
  
  #beta <- res_tmle_const$beta_unconstr
  #vcov <- res_tmle_const$cov
  #Amat <- res_tmle_const$Amat
  #w <- 1/(res_tmle_const$se^2)
  #grid <- cbind(rep(1, 101), seq(0, 1, by=0.01))
  #M=5000
  
  MC_samples <- mvrnorm(M, mu = beta, Sigma = vcov)
  
  proj_MC_samples <- apply(MC_samples, MARGIN=1, FUN=function(x){
    as.numeric(coneproj::coneA(x, amat = Amat, w = w)$thetahat)
  })
  
  proj_MC_samples <- t(proj_MC_samples)
  
  logRR_fits <- apply(proj_MC_samples, MARGIN=1, FUN=function(x){
    
    #x=proj_MC_samples[1,]
    grid %*% x
    
  })
  
  CIs <- t(apply(logRR_fits, MARGIN=1, FUN=function(x){
    
    lci = quantile(x, 0.025)
    uci = quantile(x, 0.975)
    
    c(lci, uci)
    
  }))
  
  return(CIs)
  
}
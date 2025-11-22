################
# Isotonic time-varying VE calibration
# Dependencies: coneproj
#library(coneproj)

isotone_f <- function(beta_hat, vcov, grid, indices_to_monotonize, Amat){
  
  # compute estimate of f (unconstrained)
  f_hat = grid %*% beta_hat
  
  # compute wald CI (unconstrained)
  ses_f <- sqrt(diag(grid %*% vcov %*% t(grid)))
  f_lci <- f_hat - qnorm(0.975) * ses_f
  f_uci <- f_hat + qnorm(0.975) * ses_f
  
  # compute weights for PAVA
  weights<- 1/ses_f^2
  
  # indices_to_monotonize <- which(grid[,1]==0)
  
  # Amat specifies constraint matrix for only the indices you want to monotonize
  
  if(ncol(Amat)!=length(as.numeric(f_hat)[indices_to_monotonize])){
    stop("Please provide `Amat` with ncol equal to length of monotonized indices")
  }
  
  f_mono <- c(as.numeric(f_hat)[-indices_to_monotonize], 
              coneproj::coneA(y=as.numeric(f_hat)[indices_to_monotonize], 
                            amat=Amat, 
                            w = weights[indices_to_monotonize])$thetahat)
  f_mono_lci <- c(as.numeric(f_lci)[-indices_to_monotonize],
    coneproj::coneA(y=as.numeric(f_lci)[indices_to_monotonize], 
                                amat=Amat, 
                                w = weights[indices_to_monotonize])$thetahat)
  f_mono_uci <- c(as.numeric(f_uci)[-indices_to_monotonize],
                  coneproj::coneA(y=as.numeric(f_uci)[indices_to_monotonize], 
                                amat=Amat, 
                                w = weights[indices_to_monotonize])$thetahat)
  
  return(list(
    "f_mono" = f_mono,
    "f_lci" = f_mono_lci,
    "f_uci" = f_mono_uci
  ))
}

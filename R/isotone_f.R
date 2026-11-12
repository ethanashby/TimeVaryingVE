#' isotone_f
#'
#' Post-hoc projection of time-varying VE estimates.
#' 
#' The following code projects unconstrained time-varying VE estimates obtained
#' from the semiparametric logistic regression model onto a contraint cone. This
#' is especially useful for imposing shape contraints, in our case, monotonicity
#' of the time-varying VE function.
#' 
#' @param beta_hat numeric, d-dimensional vector of coefficients associated with time-since-vaccination basis
#' @param vcov matrix, dxd covariance matrix of coefficient vector
#' @param grid matrix, lxd matrix of gridpoints that you wish to evaluate the time-varying VE function at
#' @param indices_to_monotonize numeric, vector of row indices to apply the monotone projection. Typically, will want to focus on indices after which we anticipate VE will be monotone (e.g., indices >14-days post vaccination)
#' @param Amat matrix, the constraint matrix where rows define the normal vectors of the bounding hyperplanes.
#' @return list of length three containing the monotonized estimate of f and its confidence limits
#' @examples
#'
#' @export

isotone_f <- function(beta_hat, vcov, grid, indices_to_monotonize, Amat){
  
  # compute estimate of f (unconstrained)
  f_hat = grid %*% beta_hat
  
  # compute wald CI (unconstrained)
  ses_f <- sqrt(diag(grid %*% vcov %*% t(grid)))
  f_lci <- f_hat - stats::qnorm(0.975) * ses_f
  f_uci <- f_hat + stats::qnorm(0.975) * ses_f
  
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

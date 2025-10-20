################################
# Isotonitized CI function
################################

## R code: CI for c^T beta by extremizing over ellipsoid intersect cone
library(CVXR)
library(ECOSolveR, scs)

# inputs:
# beta_hat: numeric vector length m
# Sigma: m x m positive-definite covariance matrix
# cvec: numeric vector length m (linear functional, e.g. a row of Bgrid)
# alpha: e.g. 0.05
# Dmat: (m-1) x m difference matrix so that Dmat %*% beta >= 0 enforces beta_1 <= ... <= beta_m

find_CI_ellipsoid <- function(beta_hat, Sigma, cvec, alpha, Amat){

m <- length(beta_hat)
chi2_rad <- qchisq(1 - alpha, df = m)

# compute Sigma^{-1/2} conveniently:
U <- Matrix::chol(Sigma)            # U is upper triangular, t(U) %*% U = Sigma
Sigma_half_inv <- suppressWarnings(try(solve(t(U)), silent=TRUE))  # numeric matrix = Sigma^{-1/2}

if (inherits(Sigma_half_inv, "try-error")) {
  # fallback to ginv
  Sigma_half_inv <- suppressWarnings(try(MASS::ginv(t(U)), silent=TRUE))
} else {
  
}

# helper to solve max/min
solve_extreme <- function(sign = 1) {
  beta_var <- Variable(m)
  objective <- sign * (t(cvec) %*% beta_var)   # sign = +1 for max, -1 for min
  constraints <- list( norm2(Sigma_half_inv %*% (beta_var - beta_hat)) <= sqrt(chi2_rad),
                       Amat %*% beta_var >= 0 )
  prob <- Problem(Maximize(objective), constraints = constraints)
  res <- CVXR::psolve(prob, solver = "SCS")  # or "ECOS"
  val <- as.numeric(res$value) * sign  # res$value is sign*opt; convert for min by sign=-1
  beta_opt <- res$getValue(beta_var)
  list(optval = val, beta_opt = as.numeric(beta_opt))
}

# get upper and lower
up <- solve_extreme(sign = +1)$optval
low <- solve_extreme(sign = -1)$optval
return(c(low, up))
}

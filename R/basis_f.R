##############
# Script to load useful basis functions

linear_basis <- function(x) {
  obj <- list(x = x)
  class(obj) <- "linear_basis"
  return(obj)
}

# Define predict method
predict.linear_basis <- function(object, newx, ...) {
  return(as.matrix(newx))  # just returns the input
}

# Create dimension 2 basis -- vaccine plus time-since-vaccination

psi_d2 <- function(Tvec, Vvec) {
  
  basis <- linear_basis(pmax(0, Tvec-Vvec))
  
  # Vmat: n x d
  # default: psi(T-V) = V
  
  return(list(
    mat = cbind(as.numeric(Tvec - Vvec > 0), basis$x),
    basis = basis
  ))

}

psi_bs <- function(Tvec, Vvec, df) {
  
  K = df - 3
  
  cuts <- seq(1, K) / (K + 1)
  
  #knots <- quantile(pmax(0, Tvec-Vvec)[which(Tvec > Vvec)], cuts)
  
  knots <- cuts
  
  #basis
  
  basis = bs(pmax(0, Tvec-Vvec), df=df, degree = 3, knots=knots, 
             Boundary.knots=c(-0.1, max(pmax(0, Tvec-Vvec), na.rm=TRUE)), intercept=TRUE)
  
  # Vmat: n x d
  # default: psi(T-V) = V
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      as.numeric(Tvec - Vvec > 0) * basis),
    basis = basis
  ))
}

psi_d2_early <- function(Tvec, Vvec, Vvec_early){
  
  basis <- linear_basis(pmax(0, Tvec-Vvec))
  
  # Vmat: n x d
  # default: psi(T-V) = V
  
  return(list(
    mat = cbind(as.numeric(between(Tvec, Vvec_early, Vvec)),
                as.numeric(Tvec - Vvec > 0), 
                basis$x),
    basis = basis
  ))
  
}

psi_bs_early <- function(Tvec, Vvec, Vvec_early, df, left_knot = 0, first_nonboundary = 0.30, intercept){
  
  #### Note that Vvec is vaccination date + 14 days
  #### Vvec_early is vaccination date
  
  K = df - 3
  
  cuts <- seq(first_nonboundary, K / (K + 1), length.out=K)
  
  time_post = pmax(0, Tvec-Vvec) # Time since peak VE
  
  knots <- quantile(time_post[which(Tvec > Vvec)], cuts)
  
  boundary.knots = c(left_knot, max(time_post))
  
  basis <- bs(pmax(0, Tvec-Vvec), df=df, degree=3, knots=knots, Boundary.knots = boundary.knots, intercept=intercept)
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      #as.numeric(Tvec - Vvec > 0),
      as.numeric(Tvec - Vvec > 0) * as.matrix(basis)),
    basis = basis
  ))
  
}

psi_ns_early <- function(Tvec, Vvec, Vvec_early, df, left_knot = 0, first_nonboundary = 0.30, intercept){
  
  #### Note that Vvec is vaccination date + 14 days
  #### Vvec_early is vaccination date
  
  K = df - 3
  
  cuts <- seq(first_nonboundary, K / (K + 1), length.out=K)
  
  time_post = pmax(0, Tvec-Vvec) # Time since peak VE
  
  knots <- quantile(time_post[which(Tvec > Vvec)], cuts)
  
  boundary.knots = c(left_knot, max(time_post))
  
  basis <- ns(pmax(0, Tvec-Vvec), df=df, knots=knots, Boundary.knots = boundary.knots, intercept=intercept)
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      #as.numeric(Tvec - Vvec > 0),
      as.numeric(Tvec - Vvec > 0) * as.matrix(basis)),
    basis = basis
  ))
  
}


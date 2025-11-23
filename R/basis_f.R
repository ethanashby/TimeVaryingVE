#' linear_basis
#' 
#' Helper function to create linear basis
#'
#' @param x Vector of time points
#' @return An object of type 'linear_basis'
#' @export
linear_basis <- function(x) {
  obj <- list(x = x)
  class(obj) <- "linear_basis"
  return(obj)
}

#' predict.linear.basis
#' 
#' Predict method for linear basis
#'
#' @param object Object of type 'linear_basis'
#' @param newx Vector of time points you want predictions
#' @param ... Other arguments to the predict method (shouldn't need to be used)
#' @return Matrix of predicted effectiveness
#' @export
predict.linear_basis <- function(object, newx, ...) {
  return(as.matrix(newx))  # just returns the input
}

#' psi_d2
#' 
#' Create two dimensional linear basis
#'
#' @param Tvec Vector of infection dates
#' @param Vvec Vector of vaccination dates
#' @return A list containing (i) a basis matrix evaluated at the inputs Tvec and Vvec, and (ii) the 'linear_basis' object for which additional predictions can be obtained
#' @export

psi_d2 <- function(Tvec, Vvec) {
  
  basis <- linear_basis(pmax(0, Tvec-Vvec))
  
  # Vmat: n x d
  # default: psi(T-V) = V
  
  return(list(
    mat = cbind(as.numeric(Tvec - Vvec > 0), basis$x),
    basis = basis
  ))

}

#' psi_bs
#' 
#' Create cubic b-spline basis
#'
#' @param Tvec Vector of infection dates
#' @param Vvec Vector of vaccination dates
#' @param df Degrees of freedom of splines
#' @param left_boundary_knot Location of left boundary knot (to control left boundary behavior)
#' @param leftmost_interior_knot Location of leftmost interior knot
#' @param intercept Logical, whether or not to include intercept in your basis
#' @return A list containing (i) a basis matrix evaluated at the inputs Tvec and Vvec, and (ii) the 'linear_basis' object for which additional predictions can be obtained
#' @export

psi_bs <- function(Tvec, Vvec, df, left_boundary_knot, leftmost_interior_knot, intercept) {
  
  K = df - 3
  
  cuts <- seq(leftmost_interior_knot, K / (K + 1), length.out=K)
  
  time_post = pmax(0, Tvec-Vvec) # Time since peak VE
  
  knots <- stats::quantile(time_post[which(Tvec > Vvec)], cuts)
  
  boundary.knots = c(left_boundary_knot, max(time_post))
  
  basis <- splines::bs(pmax(0, Tvec-Vvec), df=df, degree=3, knots=knots, Boundary.knots = boundary.knots, intercept=intercept)
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      as.numeric(Tvec - Vvec > 0) * basis),
    basis = basis
  ))
}

#' psi_d2_early
#' 
#' Create three dimensional linear basis (with additional scalar dimension for "early" vaccination)
#'
#' @param Tvec Vector of infection dates
#' @param Vvec Vector of dates at which vaccine is presumed to achieve peak effectiveness (often vaccination date + 14 days)
#' @param Vvec_early Vector of *actual* vaccination dates
#' @return A list containing (i) a basis matrix evaluated at the inputs Tvec and Vvec, and (ii) the 'linear_basis' object for which additional predictions can be obtained
#' @export

psi_d2_early <- function(Tvec, Vvec, Vvec_early){
  
  basis <- linear_basis(pmax(0, Tvec-Vvec))
  
  # Vmat: n x d
  # default: psi(T-V) = V
  
  return(list(
    mat = cbind(as.numeric(dplyr::between(Tvec, Vvec_early, Vvec)),
                as.numeric(Tvec - Vvec > 0), 
                basis$x),
    basis = basis
  ))
  
}

#' psi_bs_early
#' 
#' Create cubic b-spline basis (with additional scalar dimension for "early" vaccination)
#'
#' @param Tvec Vector of infection dates
#' @param Vvec Vector of dates at which vaccine is presumed to achieve peak effectiveness (often vaccination date + 14 days)
#' @param Vvec_early Vector of *actual* vaccination dates
#' @param df Degrees of freedom of splines
#' @param left_boundary_knot Location of left boundary knot (to control left boundary behavior)
#' @param leftmost_interior_knot Location of leftmost interior knot
#' @param intercept Logical, whether or not to include intercept in your basis
#' @return A list containing (i) a basis matrix evaluated at the inputs Tvec and Vvec, and (ii) the 'linear_basis' object for which additional predictions can be obtained
#' @export

psi_bs_early <- function(Tvec, Vvec, Vvec_early, df, left_boundary_knot, leftmost_interior_knot, intercept){
  
  #### Note that Vvec is vaccination date + 14 days
  #### Vvec_early is vaccination date
  
  K = df - 3
  
  cuts <- seq(leftmost_interior_knot, K / (K + 1), length.out=K)
  
  time_post = pmax(0, Tvec-Vvec) # Time since peak VE
  
  knots <- stats::quantile(time_post[which(Tvec > Vvec)], cuts)
  
  boundary.knots = c(left_boundary_knot, max(time_post))
  
  basis <- splines::bs(pmax(0, Tvec-Vvec), df=df, degree=3, knots=knots, Boundary.knots = boundary.knots, intercept=intercept)
  
  return(list(
    mat = cbind(
      as.numeric(Tvec - Vvec_early > 0 & Tvec - Vvec <= 0), 
      #as.numeric(Tvec - Vvec > 0),
      as.numeric(Tvec - Vvec > 0) * as.matrix(basis)),
    basis = basis
  ))
  
}


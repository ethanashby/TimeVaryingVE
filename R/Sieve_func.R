#' sieve_partially_linear_logistic
#'
#' Function used to estimate time-varying VE using the partially linear logistic regression model adjusting for vaccine-irrelevant infections via the sieve method.
#' 
#' @param dat data.frame, containing columns J.name, V.name, T.name, and V.early.name (optional).
#' @param J.name character, column name associated with outcome type (J=0 vaccine-irrelevant, J=1 vaccine-preventable)
#' @param V.name character, column name associated with vaccination date or peak immunization date (if V.early.name is specified)
#' @param V.early.name character, OPTIONAL column name associated with "early" effectiveness (often the vaccination date itself)
#' @param T.name character, column name associated with infection date
#' @param psi_delta function, used to build the basis for time-varying VE
#' @param verbose logical, whether to print out intermediate info
#' @param ... extra arguments to use for basis creation functions
#' @return list containing parameter estimates associated with VE, SEs, vcov matrix, the varying intercept alpha(t), and the basis
#' @export

sieve_partially_linear_logistic <- function(dat, 
                                            J.name="J", 
                                            V.name="V", 
                                            V.early.name=NULL, 
                                            T.name="T", 
                                            psi_delta, 
                                            verbose=TRUE, ...){
  
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
  data=to_regress, family=stats::binomial)
  
  beta <- coef(fit)[1:d]
  se_beta <- summary(fit)$coef[1:d,2]
  cov <- vcov(fit)[1:d, 1:d]
  alpha <- coef(fit)[(d+1):length(coef(fit))]
  
  beta_mono=beta
  
  return(
    list(
      beta = beta_mono,
      beta_unconstr = beta,
      se = se_beta,
      cov = cov,
      #Amat = Amat,
      basis = base$basis,
      alpha = alpha
    )
  )
  
}

#' create_data
#'
#' Function used to generate survival data consistent with RCTs/observational studies measuring vaccine-preventable and irrelevant infections when identification assumptions hold.
#' 
#' The following code generates data under two scenarios (i) a constant VE setting and (ii) a linearly waning VE setting on the log hazard ratio scale.
#' Survival outcomes (for both vaccine-preventable and vaccine-irrelevant pathogens) are generated from hazard frailty models.
#' 
#' @param n numeric, number of participants to simulate data for
#' @param p_boost numeric, probability of vaccination at some point during study
#' @param boost_assignment character, "random" if vaccination date is not associated with risk, "vulnerable" if to prioritize vaccinating at-risk persons early
#' @param post_boost_disinhibition binary, 0 if no increase in risk post-boost, 1 if increase in risk post-boost
#' @param lambda_0_N numeric, time-averaged baseline hazard of vaccine-irrelevant infection
#' @param lambda_0_Y numeric, time-averaged baseline hazard of vaccine-preventable infection
#' @param init_VE numeric, log(HR) of infection just after vaccination
#' @param VE_wane numeric, how log(HR) varies in time-since-vaccination
#' @param Tmax numeric, duration of study
#' @param sd_frail numeric, level of risk heterogeneity, defined as standard deviation of the lognormal frailty variable
#' @return list containing elements "covars" with summary info (including id, frailties, and vaccination status), "dat_N" containing vaccine-irrelevant infetcion data, "dat_Y_const" containing vaccine-preventable infection data with constant VE, and "dat_Y_wane" containing vaccine-preventable infection data with waning VE. 
#' @examples
#' create_data(n=100, p_boost=0.70, boost_assignment="vulnerable", post_boost_disinhibition=1, lambda_0_N=0.06, lambda_0_Y=0.08, init_VE=1.3, VE_wane = -0.5, Tmax=1, sd_frail=1.3)
#' @export

create_data<-function(n, 
                      p_boost, 
                      boost_assignment = c("random", "vulnerable"), 
                      post_boost_disinhibition = c(0,1),
                      lambda_0_N, 
                      lambda_0_Y, 
                      init_VE, 
                      VE_wane, 
                      Tmax, 
                      sd_frail){
  
  #boost_assignment=match.arg(boost_assignment)
  #post_boost_disinhibition=match.arg(boost_assignment)
  
  if(boost_assignment=="random"){
    
    if(post_boost_disinhibition==0){
      
      rho=1 #fixed frailty
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      #Random boosting so leave row and column 3 as zero
      
      myCop <- copula::normalCopula(param=copula::P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- copula::mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=0, sd=sd_frail), list(min=0, max=Tmax/p_boost)), 
                    marginsIdentical=FALSE)
      
      vars <- copula::rMvdc(n=n, myMvd)
      
    }else if(post_boost_disinhibition==1){
      
      rho=0.2 #weakly correlated change in frailties post-boost
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      #Random boosting so leave row and column 3 as zero
      
      myCop <- copula::normalCopula(param=copula::P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- copula::mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=1, sd=sd_frail), list(min=0, max=Tmax/p_boost)), 
                    marginsIdentical=FALSE)
      
      vars <- copula::rMvdc(n=n, myMvd)
      
    }
    
  }else if(boost_assignment=="vulnerable"){
    
    if(post_boost_disinhibition==0){
      
      rho=1
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      # Strong (inverse) correlation between risk and boost time: highest risk vaccinated soonest
      cor_mat[1,3]<-cor_mat[3,2]<-cor_mat[2,3]<-cor_mat[3,1]<- -0.7
      
      myCop <- copula::normalCopula(param=copula::P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- copula::mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), 
                                       list(mean=0, sd=sd_frail), 
                                       list(min=0, max=Tmax/p_boost)), 
                    marginsIdentical=FALSE)
      
      vars <- copula::rMvdc(n=n, myMvd)
      
    }else if(post_boost_disinhibition==1){
      
      rho=0.2
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      # Strong (inverse) correlation between risk and boost time: highest risk vaccinated soonest
      cor_mat[1,3]<-cor_mat[3,2]<-cor_mat[2,3]<-cor_mat[3,1]<--0.7
      
      myCop <- copula::normalCopula(param=copula::P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- copula::mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=1, sd=sd_frail), list(min=0, max=Tmax/p_boost)), 
                    marginsIdentical=FALSE)
      
      vars <- copula::rMvdc(n=n, myMvd)
      
    }
    
  }
  
  covars<-data.frame(
    id=1:n,
    U0 = vars[,1],
    U1 = vars[,2],
    Z = ifelse(vars[,3] <= Tmax, vars[,3], Inf)
  )
  
  # set hazard parameters
  
  betas_N <- c("beta0"=1, "theta0"=0, "theta1"=0, 'lambda_0'=lambda_0_N)
  betas_const_Y <- c("beta0"=1, "theta0"= init_VE, "theta1"=0, 'lambda_0'=lambda_0_Y)
  betas_wane_Y <- c("beta0"=1, "theta0"= init_VE, "theta1"=VE_wane, 'lambda_0'=lambda_0_Y)
  
  hazard_Y_2 = function(t, x, betas){
    lambda_0 = betas[["lambda_0"]]
    beta0= betas[["beta0"]]
    theta0= betas[["theta0"]]
    theta1= betas[["theta1"]]
    
    t_trt = x[["Z"]]
    tau = pmax(0, t - t_trt)
    U0 = x[["U0"]]
    U1 = x[["U1"]]
    
    return(lambda_0 * (1-0.5*sin(t * 2 * pi)) * exp(theta0 * I(t-t_trt >=0) + 
                                                      theta1 * I(t-t_trt >=0) * tau + 
                                                      beta0 * U0 * I(t-t_trt < 0) +
                                                      beta0 * U1 * I(t-t_trt >= 0)))
    
  }
  
  hazard_N = function(t, x, betas){
    lambda_0 = betas[["lambda_0"]]
    beta0= betas[["beta0"]]
    theta0= betas[["theta0"]]
    theta1= betas[["theta0"]]
    
    t_trt = x[["Z"]]
    tau = pmax(0, t - t_trt)
    U0 = x[["U0"]]
    U1 = x[["U1"]]
    
    return(lambda_0 * exp(beta0 * U0 * I(t-t_trt < 0) + beta0 * U1 * I(t-t_trt >= 0)))
  }
  
  # Simulate survival times
  
  dat_N <- simsurv::simsurv(betas = betas_N, 
                            x = covars,
                            hazard = hazard_N,
                            maxt = Tmax,
                            rootfun = stats::qlogis,
                            ids = covars$id,
                            idvar = "id")
  
  dat_Y_const <- simsurv::simsurv(betas = betas_const_Y, 
                                  x = covars,
                                  hazard = hazard_Y_2,
                                  maxt = Tmax,
                                  rootfun = stats::qlogis,
                                  ids = covars$id,
                                  idvar = "id")
  
  dat_Y_wane <- simsurv::simsurv(betas = betas_wane_Y, 
                                 x = covars,
                                 hazard = hazard_Y_2,
                                 maxt = Tmax,
                                 rootfun = stats::qlogis,
                                 ids = covars$id,
                                 idvar = "id")
  
  #return survival times
  return(
    list(
      "covars"=covars,
      "dat_N"=dat_N,
      "dat_Y_const"=dat_Y_const,
      "dat_Y_wane"=dat_Y_wane
    )
  )
  
}


#'
#' create_data_violations
#'
#' Function used to generate survival data for RCTs/observational studies
#' measuring vaccine-preventable (Y) and vaccine-irrelevant/negative-control
#' (N) infections when one or more identification assumptions may be violated.
#'
#' The function accommodates:
#' (i) a non-zero treatment effect on the negative-control endpoint,
#' (ii) a shift in circulating variants during follow-up, such that vaccine
#'      effectiveness is a time-varying interpolation between effectiveness
#'      against strains A and B, and
#' (iii) violation of the equiconfounding assumption through differential
#'      effects of unmeasured frailty on the two endpoints.
#'
#' @param n numeric, number of participants to simulate
#' @param p_boost numeric, probability of vaccination at some point during study
#' @param boost_assignment character, "random" or "vulnerable"
#' @param post_boost_disinhibition binary, 0 or 1
#' @param variant character, one of either "single" (for only simulating a 
#'        single variant) or "shift" (for simulating a shift from variant A to 
#'        variant B during the study)
#' @param lambda_0_N numeric, baseline hazard of vaccine-irrelevant infection
#' @param lambda_0_Y_A numeric, baseline hazard of vaccine-preventable
#'        infection due to strain A
#' @param lambda_0_Y_B numeric, baseline hazard of vaccine-preventable
#'        infection due to strain B
#' @param init_VE_A numeric, log(HR) for strain A immediately after vaccination
#' @param init_VE_B numeric, log(HR) for strain B immediately after vaccination
#' @param VE_wane_A numeric, change in log(HR) per unit time since vaccination
#'        for strain A
#' @param VE_wane_B numeric, change in log(HR) per unit time since vaccination
#'        for strain B
#' @param variant_midpoint numeric, calendar time at which the circulating
#'        variant prevalence is 50\%
#' @param variant_steepness numeric, controls how rapidly the circulating
#'        variant shifts from strain A to strain B
#' @param nco_hr numeric, hazard ratio for vaccination on the NCO endpoint;
#'        nco_hr = 1 corresponds to no treatment effect
#' @param bias_equiv_factor numeric, ratio of the effect of post-vaccination
#'        frailty on the NCO endpoint relative to the primary endpoint;
#'        1 corresponds to equiconfounding
#' @param Tmax numeric, duration of study
#' @param baseline_vaxx logical, denotes whether vaccination should occur 
#'        uniformly between 0 and Tmax (default) or occur at baseline 
#'        with probability 1
#' @param sd_frail numeric, standard deviation of the lognormal frailty
#'
#' @return list containing "covars", "dat_N", "dat_Y_const", and
#' "dat_Y_variant".
#'
#' @export
#'

create_data_violations <- function(
    n,
    p_boost,
    boost_assignment = c("random", "vulnerable"),
    post_boost_disinhibition = c(0, 1),
    variant = c("single", "shift"),
    lambda_0_N,
    lambda_0_Y_A,
    lambda_0_Y_B,
    init_VE_A,
    init_VE_B,
    VE_wane_A = 0,
    VE_wane_B = 0,
    Tmax,
    baseline_vaxx = FALSE,
    variant_midpoint = Tmax / 2,
    variant_steepness = 10,
    nco_hr = 1,
    bias_equiv_factor = 1,
    sd_frail
) {
  
  variant <- match.arg(variant)
  boost_assignment <- match.arg(boost_assignment)
  stopifnot(post_boost_disinhibition %in% c(0, 1))
  stopifnot(nco_hr > 0)
  stopifnot(bias_equiv_factor > 0)
  
  # ------------------------------------------------------------
  # 1. Generate frailties and vaccination times
  # ------------------------------------------------------------
  
  if (boost_assignment == "random") {
    
    if (post_boost_disinhibition == 0) {
      
      rho <- 1
      
      cor_mat <- matrix(0, nrow = 3, ncol = 3)
      diag(cor_mat) <- 1
      cor_mat[1, 2] <- cor_mat[2, 1] <- rho
      
      myCop <- copula::normalCopula(
        param = copula::P2p(cor_mat),
        dim = 3,
        dispstr = "un"
      )
      
      myMvd <- copula::mvdc(
        copula = myCop,
        margins = c(rep("norm", 2), "unif"),
        paramMargins = list(
          list(mean = 0, sd = sd_frail),
          list(mean = 0, sd = sd_frail),
          list(
            min = 0,
            max = Tmax/p_boost
          )
        ),
        marginsIdentical = FALSE
      )
      
      vars <- copula::rMvdc(n = n, myMvd)
      
    } else if (post_boost_disinhibition == 1) {
      
      rho <- 0.2
      
      cor_mat <- matrix(0, nrow = 3, ncol = 3)
      diag(cor_mat) <- 1
      cor_mat[1, 2] <- cor_mat[2, 1] <- rho
      
      myCop <- copula::normalCopula(
        param = copula::P2p(cor_mat),
        dim = 3,
        dispstr = "un"
      )
      
      myMvd <- copula::mvdc(
        copula = myCop,
        margins = c(rep("norm", 2), "unif"),
        paramMargins = list(
          list(mean = 0, sd = sd_frail),
          list(mean = 1, sd = sd_frail),
          list(
            min = 0,
            max = Tmax/p_boost
          )
        ),
        marginsIdentical = FALSE
      )
      
      vars <- copula::rMvdc(n = n, myMvd)
    }
    
  } else if (boost_assignment == "vulnerable") {
    
    if (post_boost_disinhibition == 0) {
      
      rho <- 1
      
      cor_mat <- matrix(0, nrow = 3, ncol = 3)
      diag(cor_mat) <- 1
      cor_mat[1, 2] <- cor_mat[2, 1] <- rho
      
      # Highest-risk individuals are vaccinated earliest
      cor_mat[1, 3] <- cor_mat[3, 2] <-
        cor_mat[2, 3] <- cor_mat[3, 1] <- -0.7
      
      myCop <- copula::normalCopula(
        param = copula::P2p(cor_mat),
        dim = 3,
        dispstr = "un"
      )
      
      myMvd <- copula::mvdc(
        copula = myCop,
        margins = c(rep("norm", 2), "unif"),
        paramMargins = list(
          list(mean = 0, sd = sd_frail),
          list(mean = 0, sd = sd_frail),
          list(
            min = 0,
            max = Tmax/p_boost
          )
        ),
        marginsIdentical = FALSE
      )
      
      vars <- copula::rMvdc(n = n, myMvd)
      
    } else if (post_boost_disinhibition == 1) {
      
      rho <- 0.2
      
      cor_mat <- matrix(0, nrow = 3, ncol = 3)
      diag(cor_mat) <- 1
      cor_mat[1, 2] <- cor_mat[2, 1] <- rho
      
      cor_mat[1, 3] <- cor_mat[3, 2] <-
        cor_mat[2, 3] <- cor_mat[3, 1] <- -0.7
      
      myCop <- copula::normalCopula(
        param = copula::P2p(cor_mat),
        dim = 3,
        dispstr = "un"
      )
      
      myMvd <- copula::mvdc(
        copula = myCop,
        margins = c(rep("norm", 2), "unif"),
        paramMargins = list(
          list(mean = 0, sd = sd_frail),
          list(mean = 1, sd = sd_frail),
          list(
            min = 0,
            max = Tmax/p_boost
          )
        ),
        marginsIdentical = FALSE
      )
      
      vars <- copula::rMvdc(n = n, myMvd)
    }
  }
  
  covars <- data.frame(
    id = 1:n,
    U0 = vars[, 1],
    U1 = vars[, 2],
    Z = ifelse(vars[, 3] <= Tmax, vars[, 3], Inf)
  )
  
  if(baseline_vaxx == TRUE){
  
  covars$Z <- ifelse(is.infinite(covars$Z), Inf, 0)
    
  }
  
  
  # ------------------------------------------------------------
  # 2. Hazard parameters
  # ------------------------------------------------------------
  
  betas_N <- c(
    beta0 = 1,
    theta0 = log(nco_hr),
    theta1 = 0,
    lambda_0 = lambda_0_N
  )
  
  betas_Y_A <- c(
    beta0 = 1,
    theta0 = init_VE_A,
    theta1 = VE_wane_A,
    lambda_0 = lambda_0_Y_A
  )
  
  betas_Y_B <- c(
    beta0 = 1,
    theta0 = init_VE_B,
    theta1 = VE_wane_B,
    lambda_0 = lambda_0_Y_B
  )
  
  
  # ------------------------------------------------------------
  # 3. Variant prevalence over calendar time
  # ------------------------------------------------------------
  
  p_B <- function(t, variant) {
    
    if(variant=="shift"){
    
      return(
        stats::plogis(
      variant_steepness *
        (t - variant_midpoint)
        )
    )}else{
    
      return(0)
      
    }
    
  }
  
  
  # ------------------------------------------------------------
  # 4. Primary endpoint hazard
  # ------------------------------------------------------------
  
  hazard_Y_components <- function(t, x, betas_A, betas_B) {
    
    t_trt <- x[["Z"]]
    
    tau <- pmax(0, t - t_trt)
    
    U0 <- x[["U0"]]
    U1 <- x[["U1"]]
    
    # Variant B prevalence at calendar time t
    pB <- p_B(t, variant)
    
    # Vaccine HRs
    hr_A <- exp(
      betas_A[["theta0"]] +
        betas_A[["theta1"]] * tau
    )
    
    hr_B <- exp(
      betas_B[["theta0"]] +
        betas_B[["theta1"]] * tau
    )
    
    # No vaccine effect before vaccination
    hr_A <- ifelse(t >= t_trt, hr_A, 1)
    hr_B <- ifelse(t >= t_trt, hr_B, 1)
    
    # Frailty
    frailty <- exp(
      U0 * (t < t_trt) +
        U1 * (t >= t_trt)
    )
    
    # Variant-specific hazards
    lambda_A <- frailty *
      (1-0.5*sin(t * 2 * pi)) * betas_A[["lambda_0"]] *
      hr_A
    
    lambda_B <- frailty *
      (1-0.5*sin(t * 2 * pi)) * betas_B[["lambda_0"]] *
      hr_B
    
    # Mixture contributions
    lambda_A_mixture <- (1 - pB) * lambda_A
    lambda_B_mixture <- pB * lambda_B
    
    list(
      lambda_A = lambda_A,
      lambda_B = lambda_B,
      lambda_A_mixture = lambda_A_mixture,
      lambda_B_mixture = lambda_B_mixture,
      lambda_total =
        lambda_A_mixture + lambda_B_mixture,
      p_B = pB
    )
  }
  
  hazard_Y_variant <- function(t, x, betas, variant) {
    
    components <- hazard_Y_components(
      t = t,
      x = x,
      betas_A = betas_Y_A,
      betas_B = betas_Y_B
    )
    
    components$lambda_total
  }
  
  # ------------------------------------------------------------
  # 5. Negative-control endpoint hazard
  # ------------------------------------------------------------
  
  hazard_N <- function(t, x, betas) {
    
    t_trt <- x[["Z"]]
    
    U0 <- x[["U0"]]
    U1 <- x[["U1"]]
    
    return(
      betas[["lambda_0"]] *
        exp(
          betas[["beta0"]] *
            U0 * (t < t_trt) +
            betas[["beta0"]] *
            U1 * (t >= t_trt)
        ) *
        ifelse(t >= t_trt,
               exp(betas[["theta0"]]),
               1)
    )
  }
  
  
  # ------------------------------------------------------------
  # 6. Modify NCO frailty effect to violate equiconfounding
  # ------------------------------------------------------------
  
  hazard_N_violation <- function(t, x, betas) {
    
    t_trt <- x[["Z"]]
    
    U0 <- x[["U0"]]
    U1 <- x[["U1"]]
    
    return(
      betas[["lambda_0"]] *
        exp(
          bias_equiv_factor *
            betas[["beta0"]] *
            U0 * (t < t_trt) +
            bias_equiv_factor *
            betas[["beta0"]] *
            U1 * (t >= t_trt)
        ) *
        ifelse(
          t >= t_trt,
          exp(betas[["theta0"]]),
          1
        )
    )
  }
  
  
  # ------------------------------------------------------------
  # 7. Simulate survival times
  # ------------------------------------------------------------
  
  dat_N <- simsurv::simsurv(
    betas = betas_N,
    x = covars,
    hazard = hazard_N_violation,
    maxt = Tmax,
    rootfun = stats::qlogis,
    ids = covars$id,
    idvar = "id"
  )
  
  
  dat_Y <- simsurv::simsurv(
    betas = c(
      betas_Y_A,
      betas_Y_B
    ),
    x = covars,
    hazard = hazard_Y_variant,
    maxt = Tmax,
    rootfun = stats::qlogis,
    ids = covars$id,
    idvar = "id"
  )
  
  dat_Y$variant <- NA_character_
  
  event_idx <- which(dat_Y$status == 1)
  
  if(variant=="single"){
    
    dat_Y$variant[which(dat_Y$status==1)] <- "A"
    
  }else{
    
    for (i in event_idx) {
    
    # Participant's simulated event time
    t_event <- dat_Y$eventtime[i]
    
    # Participant's covariates
    x_i <- covars[covars$id == dat_Y$id[i], ]
    
    # Evaluate the two competing hazards at the
    # participant's actual event time
    components <- hazard_Y_components(
      t = t_event,
      x = x_i,
      betas_A = betas_Y_A,
      betas_B = betas_Y_B
    )
    
    # Conditional probability that the event was caused by B
    prob_B <-
      components$lambda_B_mixture /
      components$lambda_total
    
    # Assign cause
    dat_Y$variant[i] <- ifelse(
      runif(1) < prob_B,
      "B",
      "A"
    )
    }
  }
  # ------------------------------------------------------------
  # 8. Return data
  # ------------------------------------------------------------
  
  return(
    list(
      covars = covars,
      dat_N = dat_N,
      dat_Y = dat_Y
    )
  )
}
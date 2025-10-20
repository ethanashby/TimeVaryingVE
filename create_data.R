##################################################################
# Generate test dataset to evaluate time-varying VE functions
##################################################################

# Dependencies
library(tidyverse)
library(simsurv)
library(copula)

# Function to create survival data

create_data<-function(n, p_boost, boost_assignment, post_boost_disinhibition, lambda_0_N, lambda_0_Y, init_VE, VE_wane, Tmax, sd_frail){
  
  ### Simulate covariates that determine hazard
  
  #n=1000
  #rho = 0.2
  #Tmax=1
  #boost_assignment="vulnerable"
  #post_boost_disinhibition=1
  #U1_shift = 1
  #p_boost=0.90
  #init_VE = -1
  #VE_wane = 1
  #lambda_0_Y=0.05
  #lambda_0_N=0.038
  
  if(boost_assignment=="random"){
    
    if(post_boost_disinhibition==0){
      
      rho=1 #fixed frailty
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      #Random boosting so leave row and column 3 as zero
      
      myCop <- normalCopula(param=P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=0, sd=sd_frail), list(min=0, max=Tmax + (1-p_boost))), 
                    marginsIdentical=FALSE)
      
      vars <- rMvdc(n, myMvd)
      
    }else if(post_boost_disinhibition==1){
      
      rho=0.2 #weakly correlated change in frailties post-boost
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      #Random boosting so leave row and column 3 as zero
      
      myCop <- normalCopula(param=P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=1, sd=sd_frail), list(min=0, max=Tmax + (1-p_boost))), 
                    marginsIdentical=FALSE)
      
      vars <- rMvdc(n, myMvd)
      
    }
    
  }else if(boost_assignment=="vulnerable"){
    
    if(post_boost_disinhibition==0){
      
      rho=1
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      # Strong (inverse) correlation between risk and boost time: highest risk vaccinated soonest
      cor_mat[1,3]<-cor_mat[3,2]<-cor_mat[2,3]<-cor_mat[3,1]<--0.7
      
      myCop <- normalCopula(param=P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=0, sd=sd_frail), list(min=0, max=Tmax + (1-p_boost))), 
                    marginsIdentical=FALSE)
      
      vars <- rMvdc(n, myMvd)
      
    }else if(post_boost_disinhibition==1){
      
      rho=0.2
      cor_mat = matrix(0, nrow=3, ncol=3)
      diag(cor_mat)<-1
      cor_mat[1,2]<-cor_mat[2,1]<-rho
      # Strong (inverse) correlation between risk and boost time: highest risk vaccinated soonest
      cor_mat[1,3]<-cor_mat[3,2]<-cor_mat[2,3]<-cor_mat[3,1]<--0.7
      
      myCop <- normalCopula(param=P2p(cor_mat), dim = 3, dispstr = "un")
      
      myMvd <- mvdc(copula=myCop, margins=c(rep("norm", 2), "unif"),
                    paramMargins= list(list(mean=0, sd=sd_frail), list(mean=1, sd=sd_frail), list(min=0, max=Tmax + (1-p_boost))), 
                    marginsIdentical=FALSE)
      
      vars <- rMvdc(n, myMvd)
      
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
                            rootfun = qlogis,
                            ids = covars$id,
                            idvar = "id")
  
  dat_Y_const <- simsurv::simsurv(betas = betas_const_Y, 
                                  x = covars,
                                  hazard = hazard_Y_2,
                                  maxt = Tmax,
                                  rootfun = qlogis,
                                  ids = covars$id,
                                  idvar = "id")
  
  dat_Y_wane <- simsurv::simsurv(betas = betas_wane_Y, 
                                 x = covars,
                                 hazard = hazard_Y_2,
                                 maxt = Tmax,
                                 rootfun = qlogis,
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

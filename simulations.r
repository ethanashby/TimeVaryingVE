library(SimEngine,
        lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
library(parallel)
library(tidyverse)
library(survival)
library(simsurv,
        lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
library(Hmisc)
library(coneproj,
        lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
library(quadprog)
library(mvtnorm)
library(progress)
library(copula)
library(mgcv)
library(splines)
library(msm, lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag

source("~/Desktop/GitHub/TimeVaryingVE/basis_f.R")
source("~/Desktop/GitHub/TimeVaryingVE/isotone_f.R")
source("~/Desktop/GitHub/TimeVaryingVE/create_data.R")
source("~/Desktop/GitHub/TimeVaryingVE/TMLE_func.R")
source("~/Desktop/GitHub/TimeVaryingVE/Sieve_func.R")

n_cores <- 800

SimEngine::run_on_cluster(
  
  first={
    library(SimEngine,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(parallel)
    library(tidyverse)
    library(survival)
    library(simsurv,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(Hmisc)
    library(coneproj,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(quadprog)
    library(mvtnorm)
    library(progress)
    library(copula)
    library(mgcv)
    library(splines)
    library(msm, lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    
    options(future.globals.maxSize = 8000 * 1024^2)
    
    sim <- SimEngine::new_sim()
    
    sim %<>% SimEngine::set_levels(
      n=c(10000), #sample size
      #p_boost = c(0.85, 1), #proportion boosted
      p_boost = c(0.85),
      boost_assignment = c("random", "vulnerable"), #is boosting status random or based on vulnerability
      post_boost_disinhibition = c(0, 1), #shift in unmeasured variables after post-boost disinhibition
      lambda_0_N = c(0.03), # time-averaged baseline hazard of off-target infection
      lambda_0_Y = c(0.06), # time-averaged baseline hazard of COVID-19
      init_VE = c(-1),
      VE_wane = c(1),
      Tmax = c(1),
      #Bmax = c(0.65, 1),
      Bmax = c(1),
      sd_frail = c(1, 2)
    )
    
    #id1<-sim$levels_grid$level_id[sim$levels_grid$p_boost==0.85 & sim$levels_grid$Bmax==1]
    #id2<-sim$levels_grid$level_id[sim$levels_grid$p_boost==1 & sim$levels_grid$Bmax==0.65]
    
    #sim %<>% set_levels(.keep=c(id1, id2))
    
    sim %<>% set_script(function() {
      
      # 1) Create data
      
      #L<-sim$levels_grid[6,]
      
      df <- create_data(n=L$n, 
                        p_boost=L$p_boost,
                        boost_assignment = L$boost_assignment,
                        post_boost_disinhibition=L$post_boost_disinhibition,
                        lambda_0_N=L$lambda_0_N, 
                        lambda_0_Y=L$lambda_0_Y, 
                        init_VE=L$init_VE, 
                        VE_wane=L$VE_wane, 
                        Tmax=L$Tmax, 
                        #Bmax=L$Bmax, 
                        sd_frail=L$sd_frail)
      
      # 2) Appropriately package data
      dat_N<-df$dat_N
      dat_Y_const<-df$dat_Y_const
      dat_Y_wane <- df$dat_Y_wane
      covars<-df$covars
      
      #######
      # 3) Analyze constant VE dataset
      #######
      
      # Full TTE dataset
      dat_Y_const<-data.frame(
        id=dat_Y_const$id,
        Y=dat_Y_const$eventtime,
        Delta_Y=dat_Y_const$status,
        V=covars$Z
      )
      
      # Case only dataset
      dat_const_run <- data.frame(
        id=dat_N$id,
        N=dat_N$eventtime,
        Delta_N=dat_N$status,
        Y=dat_Y_const$Y,
        Delta_Y=dat_Y_const$Delta_Y,
        V=covars$Z
      ) %>%
        filter(Delta_Y==1 | Delta_N==1) %>%
        transmute(
          `T` = pmin(Y, N),
          `A` = as.numeric(V <= T),
          `J` = ifelse(`T`==Y, 1, 0),
          `V` = V
        )
      
      # Grid to summarize findings
      const_summary <- data.frame(
        A = rep(1, 101),
        tau = seq(0, 1, by=0.01)
      )
      
      # A) Naive Cox regression
      # Build start-stop dataset for TVC
      
      d_long<-tmerge(
        data1 = dat_Y_const, 
        data2 = dat_Y_const,
        id = id, 
        tstop = Y,
        event = event(dat_Y_const$Y, dat_Y_const$Delta_Y)
      )
      
      jumps<-data.frame(
        id=covars$id,
        time = covars$Z
      )
      
      d_long<-tmerge(
        d_long, 
        jumps,
        id = id, 
        jump = tdc(time)
      )
      
      d_long<-left_join(d_long, covars %>% dplyr::select(id, Z), by="id")
      
      cx<-coxph(Surv(tstart, tstop, event) ~ jump + tt(V), 
                tt = function(tvacc, t, ...) {pmax(0, t - tvacc)}, data=d_long)
      
      # Remove large components
      cx$y <- NULL         # Response variable (redundant if you have original data)
      cx$linear.predictors <- NULL  # Only needed for predictions
      cx$residuals <- NULL  # Drop if not using residual diagnostics
      cx$weights <- NULL    # Remove case weights if not needed
      cx$model <- NULL      # Original data used for fitting
      cx$na.action <- NULL  # Drop NA action info
      # If not using robust variance, remove x and means
      cx$x <- NULL
      cx$means <- NULL
      # Further compression
      cx$call <- NULL  # If you don’t need to reconstruct the call
      cx$terms <- NULL  # If you’re not doing post-hoc modifications
      
      const_summary$logRR_cox_est <- as.matrix(const_summary) %*% coef(cx)
      ses_cox <- apply(as.matrix(const_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(cx) %*% x))})
      const_summary$logRR_cox_lci <- const_summary$logRR_cox_est - qnorm(0.975) * ses_cox
      const_summary$logRR_cox_uci <- const_summary$logRR_cox_est + qnorm(0.975) * ses_cox
      
      # B) Sieve -- Linear Estimator
      
      res_sieve_const<-sieve_partially_linear_logistic(dat=dat_const_run, 
                                                       psi_delta = psi_d2, 
                                                       V.early.name=NULL,
                                                       monotone=FALSE, 
                                                       verbose=FALSE)
      
      Amat = matrix(0, nrow = nrow(const_summary)-1, ncol= nrow(const_summary))
      
      for(i in 1:nrow(Amat)){
        
        Amat[i,i] <- -1
        Amat[i,(i+1)] <- 1
        
      }
      
      f_sieve <- isotone_f(res_sieve_const$beta, 
                           vcov = res_sieve_const$cov[1:length(res_sieve_const$beta), 1:length(res_sieve_const$beta)], 
                           grid = as.matrix(const_summary[,1:2]), 
                           indices_to_monotonize = 1:nrow(const_summary),
                           Amat = Amat)
      
      const_summary$logRR_sieve_est <- f_sieve$f_mono
      const_summary$logRR_sieve_lci <- f_sieve$f_lci
      const_summary$logRR_sieve_uci <- f_sieve$f_uci
      
      # C) TMLE -- Linear Estimator
      
      res_tmle_const<-tmle_iterative(dat=dat_const_run, psi_delta = psi_d2, V.early.name=NULL, monotone=FALSE, verbose=FALSE)
      
      f_tmle <- isotone_f(res_tmle_const$beta, 
                           vcov = res_tmle_const$cov[1:length(res_tmle_const$beta), 1:length(res_tmle_const$beta)], 
                           grid = as.matrix(const_summary[,1:2]), 
                           indices_to_monotonize = 1:nrow(const_summary),
                           Amat = Amat)
      
      const_summary$logRR_tmle_est <- f_tmle$f_mono
      const_summary$logRR_tmle_lci <- f_tmle$f_lci
      const_summary$logRR_tmle_uci <- f_tmle$f_uci
      
      
      # C) TMLE -- Linear Estimator
      
      const_results <- const_summary %>%
        pivot_longer(cols=ends_with('est'), names_to = "method", names_prefix = "logRR_", values_to="VE") %>%
        mutate(method = gsub("_est", "", method)) %>% 
        pivot_longer(cols=ends_with("uci"), names_to="tmp", names_prefix = "logRR_", values_to="UCI") %>%
        pivot_longer(cols=ends_with("lci"), names_to="tmp1", names_prefix = "logRR_", values_to="LCI") %>%
        mutate(tmp = gsub("_uci", "", tmp), tmp1 = gsub("_lci", "", tmp1)) %>% 
        filter(method==tmp & tmp==tmp1) %>%
        group_by(method) %>%
        summarise(
          MSE = mean(c(VE - L$init_VE)^2),
          coverage = mean(between(L$init_VE, LCI, UCI))
        )
      
      gc()
      
      #######
      # 4) Analyze waning VE dataset
      #######
      
      # Full TTE dataset
      dat_Y_wane<-data.frame(
        id=dat_Y_wane$id,
        Y=dat_Y_wane$eventtime,
        Delta_Y=dat_Y_wane$status,
        V=covars$Z
      )
      
      # Case only dataset
      dat_wane_run <- data.frame(
        id=dat_N$id,
        N=dat_N$eventtime,
        Delta_N=dat_N$status,
        Y=dat_Y_wane$Y,
        Delta_Y=dat_Y_wane$Delta_Y,
        V=covars$Z
      ) %>%
        filter(Delta_Y==1 | Delta_N==1) %>%
        transmute(
          `T` = pmin(Y, N),
          `A` = as.numeric(V <= T),
          `J` = ifelse(`T`==Y, 1, 0),
          `V` = V
        )
      
      # Grid to summarize findings
      wane_summary <- data.frame(
        A = rep(1, 101),
        tau = seq(0, 1, by=0.01)
      )
      
      # A) Naive Cox regression
      # Build start-stop dataset for TVC
      
      d_long<-tmerge(
        data1 = dat_Y_wane, 
        data2 = dat_Y_wane,
        id = id, 
        tstop = Y,
        event = event(dat_Y_wane$Y, dat_Y_wane$Delta_Y)
      )
      
      jumps<-data.frame(
        id=covars$id,
        time = covars$Z
      )
      
      d_long<-tmerge(
        d_long, 
        jumps,
        id = id, 
        jump = tdc(time)
      )
      
      d_long<-left_join(d_long, covars %>% dplyr::select(id, Z), by="id")
      
      cx<-coxph(Surv(tstart, tstop, event) ~ jump + tt(V), 
                tt = function(tvacc, t, ...) {pmax(0, t - tvacc)}, data=d_long)
      
      # Remove large components
      cx$y <- NULL         # Response variable (redundant if you have original data)
      cx$linear.predictors <- NULL  # Only needed for predictions
      cx$residuals <- NULL  # Drop if not using residual diagnostics
      cx$weights <- NULL    # Remove case weights if not needed
      cx$model <- NULL      # Original data used for fitting
      cx$na.action <- NULL  # Drop NA action info
      # If not using robust variance, remove x and means
      cx$x <- NULL
      cx$means <- NULL
      # Further compression
      cx$call <- NULL  # If you don’t need to reconstruct the call
      cx$terms <- NULL  # If you’re not doing post-hoc modifications
      
      wane_summary$logRR_cox_est <- as.matrix(wane_summary) %*% coef(cx)
      ses_cox <- apply(as.matrix(wane_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(cx) %*% x))})
      wane_summary$logRR_cox_lci <- wane_summary$logRR_cox_est - qnorm(0.975) * ses_cox
      wane_summary$logRR_cox_uci <- wane_summary$logRR_cox_est + qnorm(0.975) * ses_cox
      
      # B) Sieve -- Linear Estimator
      
      res_sieve_wane<-sieve_partially_linear_logistic(dat=dat_wane_run, 
                                                       psi_delta = psi_d2, 
                                                       V.early.name=NULL,
                                                       monotone=FALSE, 
                                                       verbose=FALSE)
      
      Amat = matrix(0, nrow = nrow(const_summary)-1, ncol= nrow(const_summary))
      
      for(i in 1:nrow(Amat)){
        
        Amat[i,i] <- -1
        Amat[i,(i+1)] <- 1
        
      }
      
      f_sieve <- isotone_f(res_sieve_wane$beta, 
                           vcov = res_sieve_wane$cov[1:length(res_sieve_wane$beta), 1:length(res_sieve_wane$beta)], 
                           grid = as.matrix(wane_summary[,1:2]), 
                           indices_to_monotonize = 1:nrow(wane_summary),
                           Amat = Amat)
      
      wane_summary$logRR_sieve_est <- f_sieve$f_mono
      wane_summary$logRR_sieve_lci <- f_sieve$f_lci
      wane_summary$logRR_sieve_uci <- f_sieve$f_uci
      
      # C) TMLE -- Linear Estimator
      
      res_tmle_wane<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_d2, V.early.name=NULL, monotone=FALSE, verbose=FALSE)
      
      f_tmle <- isotone_f(res_tmle_wane$beta, 
                          vcov = res_tmle_wane$cov[1:length(res_tmle_wane$beta), 1:length(res_tmle_wane$beta)], 
                          grid = as.matrix(wane_summary[,1:2]), 
                          indices_to_monotonize = 1:nrow(wane_summary),
                          Amat = Amat)
      
      wane_summary$logRR_tmle_est <- f_tmle$f_mono
      wane_summary$logRR_tmle_lci <- f_tmle$f_lci
      wane_summary$logRR_tmle_uci <- f_tmle$f_uci
      
      
      # C) TMLE -- Linear Estimator
      
      wane_results <- wane_summary %>%
        pivot_longer(cols=ends_with('est'), names_to = "method", names_prefix = "logRR_", values_to="VE") %>%
        mutate(method = gsub("_est", "", method)) %>% 
        pivot_longer(cols=ends_with("uci"), names_to="tmp", names_prefix = "logRR_", values_to="UCI") %>%
        pivot_longer(cols=ends_with("lci"), names_to="tmp1", names_prefix = "logRR_", values_to="LCI") %>%
        mutate(tmp = gsub("_uci", "", tmp), tmp1 = gsub("_lci", "", tmp1)) %>% 
        filter(method==tmp & tmp==tmp1) %>%
        group_by(method) %>%
        summarise(
          MSE = mean(c(VE - (L$init_VE + tau * L$VE_wane))^2),
          coverage = mean(between((L$init_VE + tau * L$VE_wane), LCI, UCI))
        )
      
      gc()
      
      return (list(
        ".complex" = list(
          'const_summary'= const_summary,
          'wane_summary'= wane_summary,
          'const_results' = const_results,
          'wane_results' = wane_results
        ))
      )
      
    })
    
    sim %<>% set_config(
      num_sim = 1000,
      seed=47,
      stop_at_error=FALSE,
      n_cores = n_cores,
      packages = c("copula", "survival", "progress", "simsurv", "coneproj", "splines", "Hmisc", "quadprog", "mvtnorm", "tidyverse", "msm")
    )
  },
  
  main={
    
    source("/home/bzhang3/Ashby_Ethan/VEDebias/VEWaneF_giz.R")
    
    library(SimEngine,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(parallel)
    library(tidyverse)
    library(survival)
    library(simsurv,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(Hmisc)
    library(coneproj,
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    library(quadprog)
    library(mvtnorm)
    library(progress)
    library(copula)
    library(msm, 
            lib.loc=c("/home/bzhang3/Ashby_Ethan/R_lib/")) #flag
    
    sim %<>% SimEngine::run()
  },
  
  last = {
    #save(sim, file="/home/bzhang3/Ashby_Ethan/VEDebias/sim_res.rds")
  },
  
  cluster_config = list(js="slurm")
)



L<-list(n=16000, p_boost=0.80, boost_assignment="vulnerable", post_boost_disinhibition=1, lambda_0_N=0.05, lambda_0_Y=0.08, init_VE=-1, VE_wane=1, Tmax=1)

# Create dataset
set.seed(4747)

df <- create_data(n=L$n, 
                  p_boost=L$p_boost,
                  boost_assignment = L$boost_assignment,
                  post_boost_disinhibition=L$post_boost_disinhibition,
                  lambda_0_N=L$lambda_0_N, lambda_0_Y=L$lambda_0_Y, init_VE=L$init_VE, VE_wane=L$VE_wane, Tmax=L$Tmax)

### Package data appropriately

dat_N<-df$dat_N
dat_Y_const<-df$dat_Y_const
dat_Y_wane <- df$dat_Y_wane
covars<-df$covars

#########
#### Analyze constant VE dataset
#########

dat_const<-data.frame(
  id=dat_N$id,
  N=dat_N$eventtime,
  Delta_N=dat_N$status,
  Y=dat_Y_const$eventtime,
  Delta_Y=dat_Y_const$status,
  Z=covars$Z
)

# 1) Naive Cox

dat_const <- dat_const %>% mutate(
  A = as.numeric(Z <= Y),
  tau = pmax(0, Y - Z)
)

fit_cox_const<-coxph(Surv(Y, Delta_Y) ~ A + tt(Z),
                     data=dat_const,
                     tt=function(x, t, ...){I(x <= t) * pmax(0, t-x)})

#coef(fit_cox_const)

# 2) TMLE -- Linear Estimator

dat_const_run <- dat_const %>%
  filter(Delta_Y==1 | Delta_N==1) %>%
  transmute(
    `T` = pmin(Y, N),
    `A` = as.numeric(Z <= T),
    `J` = ifelse(`T`==Y, 1, 0),
    `V` = Z
  )

const_summary <- data.frame(
  A = rep(1, 101),
  tau = seq(0, 1, by=0.01)
)


# Cox summary

const_summary$logRR_cox_est <- as.matrix(const_summary) %*% coef(fit_cox_const)
ses_cox <- apply(as.matrix(const_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(fit_cox_const) %*% x))})
const_summary$logRR_cox_lci <- const_summary$logRR_cox_est - qnorm(0.975) * ses_cox
const_summary$logRR_cox_uci <- const_summary$logRR_cox_est + qnorm(0.975) * ses_cox

# Sieve summary

res_sieve_const<-sieve_partially_linear_logistic(dat=dat_const_run, psi_delta = psi_d2, monotone=TRUE)

const_summary$logRR_sieve_est <- as.matrix(const_summary[,1:2]) %*% res_sieve_const$beta

CI <- monotone_CI_MC(beta=res_sieve_const$beta_unconstr, 
                     vcov=res_sieve_const$cov, 
                     Amat=res_sieve_const$Amat, 
                     w = 1/(res_sieve_const$se)^2, 
                     grid = as.matrix(const_summary[,1:2]),
                     M=10000)

const_summary$logRR_sieve_lci <- CI[,1]
const_summary$logRR_sieve_uci <- CI[,2]
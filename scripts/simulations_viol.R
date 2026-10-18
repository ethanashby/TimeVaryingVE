home_path = "/home/eashby/"

library(SimEngine,
        lib.loc=paste0(home_path, "R_lib/")) #flag
library(parallel)
library(tidyverse)
library(survival)
library(simsurv,
        lib.loc=paste0(home_path, "R_lib/")) #flag
library(Hmisc)
library(coneproj,
        lib.loc=paste0(home_path, "R_lib/")) #flag
library(quadprog)
library(mvtnorm)
library(progress)
library(copula)
library(mgcv)
library(splines)
library(splines2, lib.loc=paste0(home_path, "R_lib/"))
library(msm, lib.loc=paste0(home_path, "R_lib/")) #flag

source(paste0(home_path, "VEDebias/basis_f.R"))
source(paste0(home_path, "VEDebias/isotone_f.R"))
source(paste0(home_path, "VEDebias/create_data.R"))
source(paste0(home_path, "VEDebias/TMLE_func.R"))
source(paste0(home_path, "VEDebias/Sieve_func.R"))

n_cores <- 600 #600

SimEngine::run_on_cluster(
  
  first={
    library(SimEngine,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(parallel)
    library(tidyverse)
    library(survival)
    library(simsurv,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(Hmisc)
    library(coneproj,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(quadprog)
    library(mvtnorm)
    library(progress)
    library(copula)
    library(mgcv)
    library(splines)
    library(splines2, lib.loc=paste0(home_path, "R_lib/"))
    library(msm, lib.loc=paste0(home_path, "R_lib/")) #flag
    
    options(future.globals.maxSize = 8000 * 1024^2)
    
    sim <- SimEngine::new_sim()
    
    sim %<>% SimEngine::set_levels(
      assumption_viol = c("A2", "A3", "A4"),
      level_viol = c(1,2,3)
    )

    sim %<>% set_script(function() {
      
      
      # 1) Create data
      
      if(L$assumption_viol=="A2"){
        
        df<-create_data_violations(
          n=10000,
          p_boost = 0.5,
          boost_assignment="random",
          post_boost_disinhibition = 0,
          variant = "single", # No need for second variant
          lambda_0_N = 0.03,
          lambda_0_Y_A = 0.06,
          lambda_0_Y_B = 0.06,
          init_VE_A = log(0.2), #80% VE
          init_VE_B = log(1), #0% VE
          VE_wane_A = 0, #no waning
          VE_wane_B = 0,
          #variant_midpoint = Tmax/2,
          variant_steepness = NA,
          nco_hr = ifelse(L$level_viol==1, 1.1, ifelse(L$level_viol==2, 0.95, 0.75)), # HERE IS THE CHANGE
          bias_equiv_factor = 1,
          Tmax = 1,
          sd_frail = 1
        ) 
        
      }else if(L$assumption_viol=="A3"){
        
        df<-create_data_violations(
          n=10000,
          p_boost = 0.5,
          boost_assignment="random",
          post_boost_disinhibition = 0,
          variant = "shift", # THIS IS THE CHANGE
          lambda_0_N = 0.03,
          lambda_0_Y_A = 0.06,
          lambda_0_Y_B = 0.06,
          init_VE_A = log(0.20), #80% VE
          init_VE_B = ifelse(L$level_viol==1, log(0.30), ifelse(L$level_viol==2, log(0.50), log(0.80))), #70%, 50%, 20% VE for strain B
          VE_wane_A = 0, # neither wanes
          VE_wane_B = 0, # neither wanes
          Tmax = 1,
          baseline_vaxx = TRUE,
          variant_midpoint = 1000,
          nco_hr = 1, 
          bias_equiv_factor = 1,
          sd_frail = 1
        )  
        
      }else{
        
        df<-create_data_violations(
          #n=10000,
          n=2500,
          p_boost = 0.5,
          boost_assignment="random",
          post_boost_disinhibition = ifelse(L$level_viol %in% c(1,2), 0, 1),
          variant = "single", # No need for second variant
          lambda_0_N = 0.15, #0.03,
          lambda_0_Y_A = 0.25, #0.06,
          lambda_0_Y_B = 0.25, #0.06,
          init_VE_A = log(0.2), #80% VE
          init_VE_B = log(1), #0% VE
          VE_wane_A = 0, #no waning
          VE_wane_B = 0,
          #variant_midpoint = Tmax/2,
          variant_steepness = NA,
          nco_hr = 1, 
          bias_equiv_factor = ifelse(L$level_viol==1, 
                                     1.1, 
                                     ifelse(L$level_viol==2, 2, 2)), # HERE IS THE CHANGE
          Tmax = 1,
          sd_frail = ifelse(L$level_viol==1, 
                            1, 
                            ifelse(L$level_viol==2, 2, 2))
        ) 
        
      }
      
      # 2) Appropriately package data
      dat_N <- df$dat_N
      dat_Y <- df$dat_Y
      covars<-df$covars
      
      #######
      # 3) Analyze VE dataset
      #######
      
      # Full TTE dataset
      dat_Y_test<-data.frame(
        id=dat_Y$id,
        Y=dat_Y$eventtime,
        Delta_Y=dat_Y$status,
        V=covars$Z
      )
      
      # Case only dataset
      dat_run <- data.frame(
        id=dat_N$id,
        N=dat_N$eventtime,
        Delta_N=dat_N$status,
        Y=dat_Y_test$Y,
        Delta_Y=dat_Y_test$Delta_Y,
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
        data = dat_Y_test, 
        data2 = dat_Y_test,
        id = id, 
        tstop = Y,
        event = event(dat_Y_test$Y, dat_Y_test$Delta_Y)
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
      
      res_sieve_const<-sieve_partially_linear_logistic(dat=dat_run, 
                                                       psi_delta = psi_d2, 
                                                       V.early.name=NULL,
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
      
      res_tmle_const<-tmle_partially_linear_logistic(dat=dat_run, psi_delta = psi_d2, V.early.name=NULL, verbose=FALSE)
      
      f_tmle <- isotone_f(res_tmle_const$beta, 
                          vcov = res_tmle_const$cov[1:length(res_tmle_const$beta), 1:length(res_tmle_const$beta)], 
                          grid = as.matrix(const_summary[,1:2]), 
                          indices_to_monotonize = 1:nrow(const_summary),
                          Amat = Amat)
      
      const_summary$logRR_tmle_est <- f_tmle$f_mono
      const_summary$logRR_tmle_lci <- f_tmle$f_lci
      const_summary$logRR_tmle_uci <- f_tmle$f_uci
      
      # Summary
      
      if(L$assumption_viol=="A2"){
        
        const_results <- const_summary %>%
          pivot_longer(cols=ends_with('est'), names_to = "method", names_prefix = "logRR_", values_to="VE") %>%
          mutate(method = gsub("_est", "", method)) %>% 
          pivot_longer(cols=ends_with("uci"), names_to="tmp", names_prefix = "logRR_", values_to="UCI") %>%
          pivot_longer(cols=ends_with("lci"), names_to="tmp1", names_prefix = "logRR_", values_to="LCI") %>%
          mutate(tmp = gsub("_uci", "", tmp), tmp1 = gsub("_lci", "", tmp1)) %>% 
          filter(method==tmp & tmp==tmp1) %>%
          group_by(method) %>%
          dplyr::summarise(
            MSE = mean(c(VE - log(0.2))^2),
            coverage = mean(log(0.2) >= LCI & log(0.2) <= UCI)
          )
        
      }else if(L$assumption_viol=="A3"){
        
        strain_b_efficacy = ifelse(L$level_viol==1, log(0.30), ifelse(L$level_viol==2, log(0.50), log(0.80)))
        
        const_results <- const_summary %>%
          pivot_longer(cols=ends_with('est'), names_to = "method", names_prefix = "logRR_", values_to="VE") %>%
          mutate(method = gsub("_est", "", method)) %>% 
          pivot_longer(cols=ends_with("uci"), names_to="tmp", names_prefix = "logRR_", values_to="UCI") %>%
          pivot_longer(cols=ends_with("lci"), names_to="tmp1", names_prefix = "logRR_", values_to="LCI") %>%
          mutate(tmp = gsub("_uci", "", tmp), tmp1 = gsub("_lci", "", tmp1)) %>% 
          filter(method==tmp & tmp==tmp1) %>%
          group_by(method) %>%
          dplyr::summarise(
            MSE = mean(c(VE - log(0.2))^2 * I(tau <= 0.50) + c(VE - strain_b_efficacy)^2 * I(tau > 0.50)),
            coverage = mean(
              (log(0.2) * I(tau <= 0.50) + I(tau > 0.5) * strain_b_efficacy) >= LCI & 
              (log(0.2) * I(tau <= 0.50) + I(tau > 0.5) * strain_b_efficacy) <= UCI)
          )
        
        
      }else{
        
        const_results <- const_summary %>%
          pivot_longer(cols=ends_with('est'), names_to = "method", names_prefix = "logRR_", values_to="VE") %>%
          mutate(method = gsub("_est", "", method)) %>% 
          pivot_longer(cols=ends_with("uci"), names_to="tmp", names_prefix = "logRR_", values_to="UCI") %>%
          pivot_longer(cols=ends_with("lci"), names_to="tmp1", names_prefix = "logRR_", values_to="LCI") %>%
          mutate(tmp = gsub("_uci", "", tmp), tmp1 = gsub("_lci", "", tmp1)) %>% 
          filter(method==tmp & tmp==tmp1) %>%
          group_by(method) %>%
          dplyr::summarise(
            MSE = mean(c(VE - log(0.2))^2),
            coverage = mean(log(0.2) >= LCI & log(0.2) <= UCI)
          )
        
      }
      
      stats_const <- data.frame(
        "num_Y" = sum(dat_Y_test$Delta_Y),
        "num_N" = sum(dat_N$status),
        "tot_inf" = sum(dat_Y_test$Delta_Y)+sum(dat_N$status),
        "num_first_inf" = nrow(dat_run)
      )
      
      gc()
      
      
      return (list(
        ".complex" = list(
          'const_summary'= const_summary,
          'const_results' = const_results,
          'stats_const' = stats_const
        ))
      )
      
    })
    
    sim %<>% set_config(
      num_sim = 1000,
      seed=47,
      stop_at_error=FALSE,
      n_cores = n_cores,
      packages = c("copula", 
                   "survival", 
                   "progress", 
                   "simsurv", 
                   "coneproj", 
                   "splines", 
                   "mgcv", 
                   "Hmisc", 
                   "quadprog", 
                   "mvtnorm", 
                   "tidyverse", 
                   "msm")
    )
  },
  
  main={
    
    library(SimEngine,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(parallel)
    library(tidyverse)
    library(survival)
    library(simsurv,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(Hmisc)
    library(coneproj,
            lib.loc=paste0(home_path, "R_lib/")) #flag
    library(quadprog)
    library(mvtnorm)
    library(progress)
    library(copula)
    library(mgcv)
    library(splines)
    library(splines2, lib.loc=paste0(home_path, "R_lib/"))
    library(msm, lib.loc=paste0(home_path, "R_lib/")) #flag
    
    source(paste0(home_path, "VEDebias/basis_f.R"))
    source(paste0(home_path, "VEDebias/isotone_f.R"))
    source(paste0(home_path, "VEDebias/create_data.R"))
    source(paste0(home_path, "VEDebias/TMLE_func.R"))
    source(paste0(home_path, "VEDebias/Sieve_func.R"))
    
    sim %<>% SimEngine::run()
  },
  
  last = {},
  
  cluster_config = list(js="slurm")
)
#############################
# Example Data Application/Vignette
#############################

# Load packages
if(!requireNamespace("mgcv", quietly=TRUE)) install.packages("mgcv")
library(tidyverse)
library(mgcv)
library(splines)
library(coneproj)
library(splines2)
library(survival)
library(ggnewscale)
library(patchwork)
source("../TimeVaryingVE/R/create_data.R")
source("../TimeVaryingVE/R/TMLE_func.R")
source("../TimeVaryingVE/R/Sieve_func.R")
source("../TimeVaryingVE/R/isotone_f.R")

# Run function under certain set of parameters

L<-list(n=10000, 
        p_boost=0.85, 
        boost_assignment="random", 
        post_boost_disinhibition=0, 
        lambda_0_N=0.08, 
        lambda_0_Y=0.16, 
        init_VE=-1, 
        VE_wane=1, 
        Tmax=1,
        sd_frail = 2.5)

# Create dataset
set.seed(4743)

df <- create_data(n=L$n, 
                  p_boost=L$p_boost,
                  boost_assignment = L$boost_assignment,
                  post_boost_disinhibition=L$post_boost_disinhibition,
                  lambda_0_N=L$lambda_0_N, 
                  lambda_0_Y=L$lambda_0_Y, 
                  init_VE=L$init_VE, 
                  VE_wane=L$VE_wane, 
                  Tmax=L$Tmax,
                  sd_frail = L$sd_frail)

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

# dat_const <- dat_const %>% mutate(
#   A = as.numeric(Z <= Y),
#   tau = pmax(0, Y - Z)
# )
# 
# fit_cox_const<-coxph(Surv(Y, Delta_Y) ~ A + tt(Z),
#       data=dat_const,
#       tt=function(x, t, ...){I(x <= t) * pmax(0, t-x)})

dat_Y_const <- dat_Y_const %>% rename(Y=eventtime, Delta_Y = status)

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

d_long<-left_join(d_long, covars %>% dplyr::select(id, Z) %>% rename(V=Z), by="id")

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

const_summary <- data.frame(
  A = rep(1, 101),
  tau = seq(0, 1, by=0.01)
)

const_summary$logRR_cox_est <- as.matrix(const_summary) %*% coef(cx)
ses_cox <- apply(as.matrix(const_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(cx) %*% x))})
const_summary$logRR_cox_lci <- const_summary$logRR_cox_est - qnorm(0.975) * ses_cox
const_summary$logRR_cox_uci <- const_summary$logRR_cox_est + qnorm(0.975) * ses_cox

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

# Sieve summary

res_sieve_const<-sieve_partially_linear_logistic(dat=dat_const_run, psi_delta = psi_d2, monotone=FALSE)

Amat <- matrix(0, nrow=nrow(const_summary)-1, ncol=nrow(const_summary))

for(i in 1:(nrow(const_summary)-1)){
  
  Amat[i,i] = -1
  Amat[i,(i+1)] = 1
}

fit<-isotone_f(res_sieve_const$beta, 
          vcov = res_sieve_const$cov, 
          grid = as.matrix(const_summary[,1:2]), 
          indices_to_monotonize = 1:nrow(const_summary), 
          Amat = Amat)

const_summary$logRR_sieve_est <- fit$f_mono
const_summary$logRR_sieve_lci <- fit$f_lci
const_summary$logRR_sieve_uci <- fit$f_uci

# TMLE summary

res_tmle_const<-tmle_iterative(dat=dat_const_run, psi_delta = psi_d2, monotone=FALSE)

fit<-isotone_f(res_tmle_const$beta, 
               vcov = res_tmle_const$cov, 
               grid = as.matrix(const_summary[,1:2]), 
               indices_to_monotonize = 1:nrow(const_summary), 
               Amat = Amat)

const_summary$logRR_tmle_est <- fit$f_mono
const_summary$logRR_tmle_lci <- fit$f_lci
const_summary$logRR_tmle_uci <- fit$f_uci

# Truth

const_summary$logRR_truth_est <- rep(-1, 101)

# Summarize results

const_summary <- const_summary %>%
  pivot_longer(cols = c("logRR_cox_est", 
                        "logRR_cox_lci",
                        "logRR_cox_uci",
                        "logRR_sieve_est", 
                        "logRR_sieve_lci",
                        "logRR_sieve_uci",
                        "logRR_tmle_est", 
                        "logRR_tmle_lci",
                        "logRR_tmle_uci",
                        "logRR_truth_est"), 
               names_prefix = "logRR_", 
               names_sep = "_",
               names_to=c("Method", "feature"),
               values_to="logRR_value") %>%
  mutate(
    VE = 1-exp(logRR_value)
  ) %>%
  pivot_wider(
    names_from = feature, values_from=c(VE, logRR_value)
  )

p1 <- const_summary %>%
  ggplot(.)+
  geom_jitter(aes(x=tau, y=-0.5, color=Event), width=0, height=0.25, data=dat_const_run %>% 
                mutate(tau = pmax(0, T-V),
                       Event = factor(J, levels=c(0, 1), labels=c("Irrelevant", "Preventable"))), alpha=0.3)+
  scale_color_manual(values = c("gold2", "maroon"))+
  ggnewscale::new_scale_color()+
  geom_line(aes(x=tau, y=VE_est, color=Method, linetype=Method), data=const_summary, linewidth=1.25)+
  geom_ribbon(aes(x=tau, ymin=VE_lci, ymax=VE_uci, fill=Method), data=const_summary, alpha=0.25)+
  scale_color_manual(values = c("red", "purple", "blue", "black"))+
  scale_fill_manual(values = c("red", "purple", "blue", "black"))+
  scale_linetype_manual(values = c(1, 1, 1, 2))+
  scale_y_continuous(labels=scales::percent)+
  coord_cartesian(ylim = c(-0.75, 1))+
  theme_minimal()+
  labs(title = "Constant VE (63%)", 
       #caption = "n=16,000\nBoosting confounded by susceptibility\nPost-boost behavioral disinhibition",
       x= "Time since vaccination", y="Vaccine Efficacy (VE; 1-HR)")

#########
#### Analyze waning VE dataset
#########

dat_wane<-data.frame(
  id=dat_N$id,
  N=dat_N$eventtime,
  Delta_N=dat_N$status,
  Y=dat_Y_wane$eventtime,
  Delta_Y=dat_Y_wane$status,
  Z=covars$Z
)

# 1) Naive Cox

dat_Y_wane <- dat_Y_wane %>% rename(Y=eventtime, Delta_Y = status)

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

d_long<-left_join(d_long, covars %>% dplyr::select(id, Z) %>% rename(V=Z), by="id")

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

wane_summary <- data.frame(
  A = rep(1, 101),
  tau = seq(0, 1, by=0.01)
)

wane_summary$logRR_cox_est <- as.matrix(wane_summary) %*% coef(cx)
ses_cox <- apply(as.matrix(wane_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(cx) %*% x))})
wane_summary$logRR_cox_lci <- wane_summary$logRR_cox_est - qnorm(0.975) * ses_cox
wane_summary$logRR_cox_uci <- wane_summary$logRR_cox_est + qnorm(0.975) * ses_cox

# 2) TMLE -- Linear Estimator

dat_wane_run <- dat_wane %>%
  filter(Delta_Y==1 | Delta_N==1) %>%
  transmute(
    `T` = pmin(Y, N),
    `A` = as.numeric(Z <= T),
    `J` = ifelse(`T`==Y, 1, 0),
    `V` = Z
  )

# Sieve summary

res_sieve_wane<-sieve_partially_linear_logistic(dat=dat_wane_run, psi_delta = psi_d2, monotone=FALSE)

Amat <- matrix(0, nrow=nrow(wane_summary)-1, ncol=nrow(wane_summary))

for(i in 1:(nrow(wane_summary)-1)){
  
  Amat[i,i] = -1
  Amat[i,(i+1)] = 1
}

fit<-isotone_f(res_sieve_wane$beta, 
               vcov = res_sieve_wane$cov, 
               grid = as.matrix(wane_summary[,1:2]), 
               indices_to_monotonize = 1:nrow(wane_summary), 
               Amat = Amat)

wane_summary$logRR_sieve_est <- fit$f_mono
wane_summary$logRR_sieve_lci <- fit$f_lci
wane_summary$logRR_sieve_uci <- fit$f_uci

# TMLE summary

res_tmle_wane<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_d2, monotone=FALSE)

fit<-isotone_f(res_tmle_wane$beta, 
               vcov = res_tmle_wane$cov, 
               grid = as.matrix(wane_summary[,1:2]), 
               indices_to_monotonize = 1:nrow(wane_summary), 
               Amat = Amat)

wane_summary$logRR_tmle_est <- fit$f_mono
wane_summary$logRR_tmle_lci <- fit$f_lci
wane_summary$logRR_tmle_uci <- fit$f_uci

# Truth

wane_summary$logRR_truth_est <- rep(-1, 101) + seq(0,1,0.01)


# Summarize results

wane_summary <- wane_summary %>%
  pivot_longer(cols = c("logRR_cox_est", 
                        "logRR_cox_lci",
                        "logRR_cox_uci",
                        "logRR_sieve_est", 
                        "logRR_sieve_lci",
                        "logRR_sieve_uci",
                        "logRR_tmle_est", 
                        "logRR_tmle_lci",
                        "logRR_tmle_uci",
                        "logRR_truth_est"), 
               names_prefix = "logRR_", 
               names_sep = "_",
               names_to=c("Method", "feature"),
               values_to="logRR_value") %>%
  mutate(
    VE = 1-exp(logRR_value)
  ) %>%
  pivot_wider(
    names_from = feature, values_from=c(VE, logRR_value)
  )

p2 <- ggplot(wane_summary)+
  geom_jitter(aes(x=tau, y=-0.5, color=Event), width=0, height=0.25, data=dat_wane_run %>% 
                mutate(tau = pmax(0, T-V),
                       Event = factor(J, levels=c(0, 1), labels=c("Irrelevant", "Preventable"))), alpha=0.3)+
  scale_color_manual(values = c("gold2", "maroon"))+
  ggnewscale::new_scale_color()+
  geom_line(aes(x=tau, y=VE_est, color=Method, linetype=Method), data=wane_summary, linewidth=1.25)+
  geom_ribbon(aes(x=tau, ymin=VE_lci, ymax=VE_uci, fill=Method), data=wane_summary, alpha=0.25)+
  scale_color_manual(values = c("red", "purple", "blue", "black"))+
  scale_fill_manual(values = c("red", "purple", "blue", "black"))+
  scale_linetype_manual(values = c(1, 1, 1, 2))+
  coord_cartesian(ylim=c(-0.75, 1))+
  scale_y_continuous(labels=scales::percent)+
  theme_minimal()+
  labs(title = "(Log) Linear Waning VE (63% -> 0%)", 
       #caption = "n=16,000\nBoosting confounded by susceptibility\nPost-boost behavioral disinhibition",
       x= "Time since vaccination", y="Vaccine Efficacy (VE; 1-HR)")

#pdf("../VE_vignette_dos.pdf", width=9, height=6)
p1+ theme(plot.caption = element_text(hjust = 0), plot.subtitle = element_text(size=20),
          axis.title = element_text(size=18),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14))+
p2 +theme(plot.caption = element_text(hjust = 0), plot.subtitle = element_text(size=20),
          axis.title = element_text(size=18),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14))+
plot_layout(guides = "collect") & theme(legend.position = "bottom")
# labs(
#   caption = "n=16,000\nBoosting confounded by susceptibility\nPost-boost behavioral disinhibition"
# )+
# theme(plot.caption = element_text(hjust = 0))
#dev.off()
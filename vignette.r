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
source("~/Desktop/GitHub/TimeVaryingVE/create_data.R")
source("~/Desktop/GitHub/TimeVaryingVE/TMLE_func.R")
source("~/Desktop/GitHub/TimeVaryingVE/Sieve_func.R")

# Run function under certain set of parameters

L<-list(n=10000, 
        p_boost=0.80, 
        boost_assignment="vulnerable", 
        post_boost_disinhibition=1, 
        lambda_0_N=0.05, 
        lambda_0_Y=0.08, 
        init_VE=-1, 
        VE_wane=1, 
        Tmax=1,
        sd_frail = 1)

# Create dataset
set.seed(4747)

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

# TMLE summary

res_tmle_const<-tmle_iterative(dat=dat_const_run, psi_delta = psi_d2, monotone=TRUE)

const_summary$logRR_TMLE_est <- as.matrix(const_summary[,1:2]) %*% res_tmle_const$beta

CI <- monotone_CI_MC(beta=res_tmle_const$beta_unconstr, 
                     vcov=res_tmle_const$cov, 
                     Amat=res_tmle_const$Amat, 
                     w = 1/(res_tmle_const$se)^2, grid = as.matrix(const_summary[,1:2]),
                     M=10000)

const_summary$logRR_TMLE_lci <- CI[,1]
const_summary$logRR_TMLE_uci <- CI[,2]

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
                        "logRR_TMLE_est", 
                        "logRR_TMLE_lci",
                        "logRR_TMLE_uci",
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

p1 <- ggplot(const_summary)+
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
  ylim(c(-0.75, 1))+
  theme_minimal()+
  labs(title = "Constant VE (63%)", 
       #caption = "n=16,000\nBoosting confounded by susceptibility\nPost-boost behavioral disinhibition",
       x= "Time since vaccination", y="Vaccine Efficacy (VE; 1-HR)")

#########
#### Analyze constant VE dataset
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

dat_wane <- dat_wane %>% mutate(
  A = as.numeric(Z <= Y),
  tau = pmax(0, Y - Z)
)

fit_cox_wane<-coxph(Surv(Y, Delta_Y) ~ A + tt(Z),
                     data=dat_wane,
                     tt=function(x, t, ...){I(x <= t) * pmax(0, t-x)})

#coef(fit_cox_const)

# 2) TMLE -- Linear Estimator

dat_wane_run <- dat_wane %>%
  filter(Delta_Y==1 | Delta_N==1) %>%
  transmute(
    `T` = pmin(Y, N),
    `A` = as.numeric(Z <= T),
    `J` = ifelse(`T`==Y, 1, 0),
    `V` = Z
  )

wane_summary <- data.frame(
  A = rep(1, 101),
  tau = seq(0, 1, by=0.01)
)


# Cox summary

wane_summary$logRR_cox_est <- as.matrix(wane_summary) %*% coef(fit_cox_wane)
ses_cox <- apply(as.matrix(wane_summary[,1:2]), MARGIN=1, FUN=function(x){sqrt(t(x) %*% (vcov(fit_cox_wane) %*% x))})
wane_summary$logRR_cox_lci <- wane_summary$logRR_cox_est - qnorm(0.975) * ses_cox
wane_summary$logRR_cox_uci <- wane_summary$logRR_cox_est + qnorm(0.975) * ses_cox

# Sieve summary

res_sieve_wane<-sieve_partially_linear_logistic(dat=dat_wane_run, psi_delta = psi_d2, monotone=TRUE)

wane_summary$logRR_sieve_est <- as.matrix(wane_summary[,1:2]) %*% res_sieve_wane$beta

CI <- monotone_CI_MC(beta=res_sieve_wane$beta_unconstr, 
                     vcov=res_sieve_wane$cov, 
                     Amat=res_sieve_wane$Amat, 
                     w = 1/(res_sieve_wane$se)^2, 
                     grid = as.matrix(wane_summary[,1:2]),
                     M=10000)

wane_summary$logRR_sieve_lci <- CI[,1]
wane_summary$logRR_sieve_uci <- CI[,2]

# TMLE summary

res_tmle_wane<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_d2, monotone=TRUE)

wane_summary$logRR_TMLE_est <- as.matrix(wane_summary[,1:2]) %*% res_tmle_wane$beta

CI <- monotone_CI_MC(beta=res_tmle_wane$beta_unconstr, 
                     vcov=res_tmle_wane$cov, 
                     Amat=res_tmle_wane$Amat, 
                     w = 1/(res_tmle_wane$se)^2, grid = as.matrix(wane_summary[,1:2]),
                     M=10000)

wane_summary$logRR_TMLE_lci <- CI[,1]
wane_summary$logRR_TMLE_uci <- CI[,2]

# Truth

wane_summary$logRR_truth_est <- seq(-1, 0, by=0.01)

# Summarize results

wane_summary <- wane_summary %>%
  pivot_longer(cols = c("logRR_cox_est", 
                        "logRR_cox_lci",
                        "logRR_cox_uci",
                        "logRR_sieve_est", 
                        "logRR_sieve_lci",
                        "logRR_sieve_uci",
                        "logRR_TMLE_est", 
                        "logRR_TMLE_lci",
                        "logRR_TMLE_uci",
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
  ylim(c(-0.75, 1))+
  theme_minimal()+
  labs(title = "(Log) Linear Waning VE (63% -> 0%)", 
       #caption = "n=16,000\nBoosting confounded by susceptibility\nPost-boost behavioral disinhibition",
       x= "Time since vaccination", y="Vaccine Efficacy (VE; 1-HR)")

pdf("~/Desktop/VE_vignette.pdf", width=9, height=6)
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
dev.off()
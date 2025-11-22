# Debiasing time-varying, hazard-based vaccine effectiveness estimates using vaccine-irrelevant infections

This repo includes functions to estimate time-varying vaccine
effectiveness while avoiding confounding and selection bias by adjusting
for negative control, vaccine-irrelevant infections.

## Can't I just use Cox regression?

Even in randomized studies, standard Cox regression does not provide
reliable evidence of waning vaccine effectiveness due to depletion of
susceptibles bias. We develop an alternative estimation approach based
on semiparametric partially linear logistic regression to avoid bias and
correctly infer waning vaccine effectiveness estimates from
observational data. You can find our paper detailing the methodology
[here](https://arxiv.org/abs/2511.15099).

![Example data application to hypothetical randomized trial of 10,000
participants with unmeasured heterogeneity in infection
risk.](vignette/VE_vignette.jpg){width="8in"}

## How do I apply this method?

You can find a detailed script `vignette.r` which will reproduce a mock
analysis shown in the paper's supplement. Below, we walk through a
minimal example to help you design your own mock analysis.

*Step 1*: load required packages and analysis functions

``` r
library(mgcv, tidyverse, survival)
library(coneproj)
source("../TimeVaryingVE/R/create_data.R")
source("../TimeVaryingVE/R/TMLE_func.R")
source("../TimeVaryingVE/R/basis_f.R")
```

*Step 2* initialize configuration parameters and create mock datasets.

-   `n`: number of participants in study
-   `p_boost`: probability of study participants being vaccinated
-   `boost_assignment`: choose between "random" or "vulnerable" to
    determine whether vaccination date is associated with underlying
    risk, or early vaccination is associated with higher risk.
-   `lambda_0_N`: time-averaged baseline hazard of vaccine-irrelevant
    (negative control) infection
-   `lambda_0_Y`: time-averaged baseline hazard of vaccine-preventable
    (primary) infection.
-   `init_VE`: log(HR) associated with vaccination immediately after
    vaccination.
-   `VE_wane`: how log(HR) varies (linearly) in time-since-vaccination
-   `Tmax`: duration of study (years)
-   `sd_frail`: the parameter governing the level of frailty
    heterogeneity, as $\log(U) \sim N(0, \text{sd}_{frail})$

``` r
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
dat_Y_wane <- df$dat_Y_wane
covars<-df$covars

dat_wane<-data.frame(
  id=dat_N$id,
  N=dat_N$eventtime,
  Delta_N=dat_N$status,
  Y=dat_Y_wane$eventtime,
  Delta_Y=dat_Y_wane$status,
  Z=covars$Z
)
```

*Step 3* Run proposed TMLE estimator on waning dataset.

``` r
dat_wane_run <- dat_wane %>%
  filter(Delta_Y==1 | Delta_N==1) %>%
  transmute(
    `T` = pmin(Y, N),
    `A` = as.numeric(Z <= T),
    `J` = ifelse(`T`==Y, 1, 0),
    `V` = Z
  )
  
res_tmle_wane<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_d2, monotone=FALSE)
```

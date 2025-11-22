

sim<-readRDS("./Desktop/sim.rds")

for(i in 1:length(sim$results_complex)){
  
  if(i==1){
  res<-bind_cols(
    sim$results[i,],
    bind_rows(sim$results_complex[[i]]$const_results %>% mutate(setting="Constant"), 
            sim$results_complex[[i]]$wane_results %>% mutate(setting="Waning")) %>%
    pivot_wider(names_from = "method", values_from = c("MSE", "coverage")))}else{
    res <- bind_rows(res, 
                     bind_cols(
                       sim$results[i,],
                       bind_rows(sim$results_complex[[i]]$const_results %>% mutate(setting="Constant"), 
                                 sim$results_complex[[i]]$wane_results %>% mutate(setting="Waning")) %>%
                         pivot_wider(names_from = "method", values_from = c("MSE", "coverage"))))
                     
  }
  
}


res_cover <- res %>% 
  pivot_longer(cols=starts_with("coverage"), names_prefix = "coverage_", names_to = "estimator", values_to="coverage") %>%
  group_by(setting, boost_assignment, post_boost_disinhibition, lambda_0_Y, estimator, sd_frail) %>%
  dplyr::summarise(cover = mean(coverage)) %>%
  ungroup() %>% 
  mutate(
    labels = case_when(
      (sd_frail == 1 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 1,
      (sd_frail == 2 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 2,
      (sd_frail == 1 & post_boost_disinhibition==1 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 3,
      (sd_frail == 1 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="vulnerable") ~ 4,
      .default = NA)
  ) %>%
  filter(!is.na(labels)) %>% 
  mutate(cover = round(cover * 100, 1)) %>%
  pivot_wider(names_from = c("setting", "estimator"), values_from = "cover") %>%
  dplyr::select(-labels) %>%
  mutate(
    post_boost_disinhibition = factor(ifelse(post_boost_disinhibition==0, "N", "Y"), levels=c("N", "Y")),
    sd_frail = ifelse(sd_frail==1, "Low ($\\sigma_U=1$)", "High ($\\sigma_U=2$)")
  ) %>%
  rename(
    `Vaccination Assignment` = boost_assignment,
    Disinhibition = post_boost_disinhibition,
    `Force of COVID-19 infection ($\\lambda_{01}$)` = lambda_0_Y,
    `Risk Heterogeneity ($\\sigma_U$)` = `sd_frail`,
    `Cox`=`Constant_cox`,
    `SLR (Sieve)`=`Constant_sieve`,
    `SLR (TMLE)`=`Constant_sieve`,
    `Cox2`=`Waning_cox`,
    `SLR2 (Sieve)`=`Waning_sieve`,
    `SLR2 (TMLE)`=`Waning_sieve`
  ) %>%
  kable(
    format="latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Empirical coverage probabilities when $\\lambda_{01} = 0.06$ and boost probability is $p = 0.85$ uniformly over $[0,1]$. $n = 10{,}000$."
  ) %>%
  add_header_above(c(" " = 4, "Constant VE" = 3, "Linear VE" = 3)) %>%
  kable_styling(
    latex_options = c("hold_position", "striped"),
    full_width = FALSE,
    font_size = 10
  )

res_MSE <- res %>% 
  pivot_longer(cols=starts_with("MSE"), names_prefix = "MSE_", names_to = "estimator", values_to="mse") %>%
  group_by(setting, boost_assignment, post_boost_disinhibition, lambda_0_Y, estimator, sd_frail) %>%
  dplyr::summarise(mse = mean(mse)) %>%
  ungroup() %>% 
  mutate(
    labels = case_when(
      (sd_frail == 1 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 1,
      (sd_frail == 2 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 2,
      (sd_frail == 1 & post_boost_disinhibition==1 & lambda_0_Y==0.06 & boost_assignment=="random") ~ 3,
      (sd_frail == 1 & post_boost_disinhibition==0 & lambda_0_Y==0.06 & boost_assignment=="vulnerable") ~ 4,
      .default = NA)
  ) %>%
  filter(!is.na(labels)) %>% 
  mutate(mse = round(mse, 3)) %>%
  pivot_wider(names_from = c("setting", "estimator"), values_from = "mse") %>%
  dplyr::select(-labels) %>%
  mutate(
    post_boost_disinhibition = factor(ifelse(post_boost_disinhibition==0, "N", "Y"), levels=c("N", "Y")),
    sd_frail = ifelse(sd_frail==1, "Low ($\\sigma_U=1$)", "High ($\\sigma_U=2$)")
  ) %>%
  rename(
    `Vaccination Assignment` = boost_assignment,
    Disinhibition = post_boost_disinhibition,
    `Force of COVID-19 infection ($\\lambda_{01}$)` = lambda_0_Y,
    `Risk Heterogeneity ($\\sigma_U$)` = `sd_frail`,
    `Cox`=`Constant_cox`,
    `SLR (Sieve)`=`Constant_sieve`,
    `SLR (TMLE)`=`Constant_sieve`,
    `Cox2`=`Waning_cox`,
    `SLR2 (Sieve)`=`Waning_sieve`,
    `SLR2 (TMLE)`=`Waning_sieve`
  ) %>%
  kable(
    format="latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Empirical coverage probabilities when $\\lambda_{01} = 0.06$ and boost probability is $p = 0.85$ uniformly over $[0,1]$. $n = 10{,}000$."
  ) %>%
  add_header_above(c(" " = 4, "Constant VE" = 3, "Linear VE" = 3)) %>%
  kable_styling(
    latex_options = c("hold_position", "striped"),
    full_width = FALSE,
    font_size = 10
  )


res %>% 
  pivot_longer(cols=starts_with("coverage"), names_prefix = "coverage_", names_to = "estimator", values_to="coverage") %>%
  group_by(setting, boost_assignment, post_boost_disinhibition, lambda_0_Y, estimator, sd_frail) %>%
  dplyr::summarise(cover = mean(coverage)) %>%
  filter(boost_assignment=="random" & post_boost_disinhibition==0) %>%
  ungroup() %>% 
  mutate(cover = round(cover * 100, 1)) %>%
  pivot_wider(names_from = c("setting", "estimator"), values_from = "cover") %>%
  arrange(sd_frail, lambda_0_Y) %>%
  mutate(
    post_boost_disinhibition = factor(ifelse(post_boost_disinhibition==0, "N", "Y"), levels=c("N", "Y")),
    sd_frail = ifelse(sd_frail==1, "Low ($\\sigma_U=1$)", "High ($\\sigma_U=2$)")
  ) %>%
  rename(
    `Vaccination Assignment` = boost_assignment,
    Disinhibition = post_boost_disinhibition,
    `Force of COVID-19 infection ($\\lambda_{01}$)` = lambda_0_Y,
    `Risk Heterogeneity ($\\sigma_U$)` = `sd_frail`,
    `Cox`=`Constant_cox`,
    `SLR (Sieve)`=`Constant_sieve`,
    `SLR (TMLE)`=`Constant_sieve`,
    `Cox2`=`Waning_cox`,
    `SLR2 (Sieve)`=`Waning_sieve`,
    `SLR2 (TMLE)`=`Waning_sieve`
  ) %>%
  kable(
    format="latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Empirical coverage probabilities when $\\lambda_{01} = 0.06$ and boost probability is $p = 0.85$ uniformly over $[0,1]$. $n = 10{,}000$."
  ) %>%
  add_header_above(c(" " = 4, "Constant VE" = 3, "Linear VE" = 3)) %>%
  kable_styling(
    latex_options = c("hold_position", "striped"),
    full_width = FALSE,
    font_size = 10
  )


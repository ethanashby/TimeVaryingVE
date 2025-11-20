library(haven)
library(survival)
library(tidyverse)
library(data.table)
library(quadprog)
library(MASS)
library(Matrix)
library(coneproj)
library(ggnewscale)
library(patchwork)
library(splines)
library(mgcv)

ref_dt = dmy('23SEP2021')
cens_dt = dmy("05APR2022")
strata_dt=dmy("15JAN2022")
omicron_dt = dmy("15DEC2021")

source("~/Desktop/GitHub/TimeVaryingVE/Sieve_func.R")
source("~/Desktop/GitHub/TimeVaryingVE/TMLE_func.R")
source("~/Desktop/GitHub/TimeVaryingVE/basis_f.R")
source("~/Desktop/GitHub/TimeVaryingVE/isotone_f.R")

#########
### Load datasets
#########

#adttea
adttea<-read_xpt("/Volumes/trials/covpn/p3001/analysis/openlabel/qdata/adttea.xpt")

#biofire data
biofire<-read.csv("/Volumes/trials/covpn/p3001/analysis/openlabel/adata/abiofire.csv")

#primary efficacy data
pe<-fread("/Volumes/trials/covpn/p3001/analysis/openlabel/adata/P3001_BlindedPhasePartABC_PrimaryEfficacyData.csv", header = TRUE)

#variant data
GISAID <- read.table("~/Desktop/AJE Paper/Biofire analysis/pango/gisaid_variants_statistics_tsv_2024_07_13_2044/gisaid_variants_statistics.tsv", sep = "\t", header=TRUE)

# create dataset `oracle` that contains key study period start and end dates

oracle<-pe %>% 
  filter(PARAMCD=="TTCVDC1") %>% 
  mutate(URM_status = ifelse((RACE == "WHITE" | is.na(RACE)) & (ETHNIC=="NOT HISPANIC OR LATINO" | is.na(ETHNIC)), 0, 1)) %>%
  mutate(START_A = dmy(AP01SDT), 
         END_A = pmin(dmy(AP01EDT), dmy(EOSDT), dmy(EFFCODT), na.rm=TRUE),
         START_OL = case_when(
           is.na(dmy(AP02SDT)) ~ NA,
           !is.na(dmy(AP02SDT)) ~ dmy(AP02SDT)),
         START_XO = case_when(
           is.na(dmy(TR02SDT)) ~ NA,
           !is.na(dmy(TR02SDT)) ~ dmy(TR02SDT)),
         START_BOOST = case_when(
           is.na(dmy(TR03SDT)) ~ NA,
           !is.na(dmy(TR03SDT)) ~ dmy(TR03SDT)),
         END_FU = case_when(
           is.na(pmin(dmy(EOSDT), dmy("05APR2022"), dmy(EFFCODT), na.rm=TRUE)) ~ NA,
           !is.na(pmin(dmy(EOSDT), dmy("05APR2022"), dmy(EFFCODT), na.rm=TRUE)) ~ pmin(dmy(EOSDT), dmy(EFFCODT), dmy("05APR2022"), na.rm=TRUE))) %>%
  mutate(
    `Fall2020_FU`=case_when(
      END_FU< dmy('27JUL2020') ~ NA, 
      START_A>dmy("30NOV2020") ~ NA, 
      END_FU>=dmy('27JUL2020') & START_A<=dmy("30NOV2020")~ (as.numeric(difftime(pmin(END_FU, dmy("30NOV2020")), pmax(dmy('27JUL2020'), START_A), units="days"))+1)/365.25),
    `Winter2020_FU`=case_when(END_FU< dmy('01DEC2020') ~ NA, START_A>dmy("28FEB2021") ~ NA, END_FU>=dmy('01DEC2020') & START_A<=dmy("28FEB2021")~ (as.numeric(difftime(pmin(END_FU, dmy("28FEB2021")), pmax(dmy('01DEC2020'), START_A), units="days")+1))/365.25),
    `Spring2021_FU`=case_when(END_FU< dmy('01MAR2021') ~ NA, START_A>dmy("31MAY2021") ~ NA, END_FU>=dmy('01MAR2021') & START_A<=dmy("31MAY2021")~ as.numeric(difftime(pmin(END_FU, dmy("31MAY2021")), pmax(dmy('01MAR2021'), START_A), units="days")+1)/365.25), 
    `Summer2021_FU`=case_when(END_FU< dmy('01JUN2021') ~ NA, START_A>dmy("31AUG2021") ~ NA, END_FU>=dmy('01JUN2021') & START_A<=dmy("31AUG2021")~ as.numeric(difftime(pmin(END_FU, dmy("31AUG2021")), pmax(dmy('01JUN2021'), START_A), units="days")+1)/365.25), 
    `Fall2021_FU`=case_when(END_FU< dmy('01SEP2021') ~ NA, START_A>dmy("30NOV2021") ~ NA, END_FU>=dmy('01SEP2021') & START_A<=dmy("30NOV2021")~ as.numeric(difftime(pmin(END_FU, dmy("30NOV2021")), pmax(dmy('01SEP2021'), START_A), units="days")+1)/365.25), 
    `Winter2021_FU`=case_when(
      END_FU < dmy('01DEC2021') ~ NA, 
      START_A > dmy("28FEB2022") ~ NA, 
      END_FU>=dmy('01DEC2021') & START_A<=dmy("28FEB2022")~ as.numeric(difftime(pmin(END_FU, dmy("28FEB2022")), pmax(dmy('01DEC2021'), START_A), units="days")+1)/365.25),
    `Spring2022_FU`=case_when(END_FU< dmy('01MAR2022') ~ NA, START_A>dmy("05APR2022") ~ NA, END_FU>=dmy('01MAR2022') & START_A<=dmy("05APR2022")~ as.numeric(difftime(pmin(END_FU, dmy("05APR2022")), pmax(dmy('01MAR2022'), START_A), units="days")+1)/365.25))

# Filter for TTE during booster phase among PPPSFL population

adttea_e <- adttea %>% 
  filter(PARAMCD=="TTCVD24" & PPPSFL=="Y" & ADT > dmy("23SEP2021")) %>%
  dplyr::select(SUBJID,AGE,SEX,RACE,ETHNIC,STRATAR,STRATARN,TRT01AN,SAFBDFL,PPPSFL,PPBDFL,PARAMCD,TR02SDT,TR03SDT,MRNA2SDT,ADT,CNSR,FSTINFDT,FPBEVTDT,LGROUP, LINEAGE, XQDT,TRT01A,TR01SDT) %>%
  filter(ymd(ADT) > ref_dt)

# Obtain GISAID data from study period and impute missing lineage groups

GISAID_boost <- GISAID %>% 
  filter(Country=="USA" & Type=="Variant") %>% 
  filter(ymd(Week.prior.to) >= dmy("23SEP2021") & ymd(Week.prior.to) <= dmy("5APR2022")+7)

set.seed(44)
for(i in 1:nrow(adttea_e)){
  if(i%%100==0){print(i)}
  if(adttea_e$CNSR[i]==1){next}else{
    if(str_length(adttea_e$LGROUP[i])>0){next}else{
      
      indx<-GISAID_boost$Week.prior.to[which.max(ymd(adttea_e$ADT[i]) < ymd(GISAID_boost$Week.prior.to))]
      tmp<-GISAID_boost %>% filter(Week.prior.to==indx)
      
      res<-sample(tmp$Value, size = 1, prob=tmp$X..per.Country.and.Week)
      adttea_e$LGROUP[i] <- case_when(
        grepl("Omicron", res) == TRUE | grepl("omicron", res)==TRUE ~ "Omicron",
        grepl("Delta", res) == TRUE | grepl("delta", res)==TRUE ~ "Delta",
        .default = "Other"
      )
    }
  }
}

###############
# create dataset 'cove' which contains all the data you will need for SARS-CoV-2 data analysis
###############

cove<-adttea_e

# collect analysis IDs: 20,404 participants

analysis_id<-unique(cove$SUBJID)

cove<-left_join(cove, oracle %>% 
                  dplyr::select(SUBJID, standardized_risk_score, END_FU), 
                by="SUBJID")

#flag SARS-CoV-2 events
cove$SARSevent = as.numeric(cove$CNSR==0 & between(cove$ADT, ref_dt, pmin(cove$END_FU, cens_dt)))

#flag time to SARS-CoV-2 event
cove$SARStte = ifelse(
  cove$SARSevent == 0,
  as.numeric(difftime(pmin(ymd(cove$END_FU), cens_dt, na.rm=TRUE), ref_dt, units="days")),
  as.numeric(difftime(pmin(ymd(cove$END_FU), cens_dt, ymd(cove$ADT), na.rm=TRUE), ref_dt, units="days"))
)

#Calculate standardized risk score and minority indicators
cove$standardized_risk_score<-ifelse(is.na(cove$standardized_risk_score), 0, cove$standardized_risk_score)
cove$MinorityInd <- ifelse(cove$RACE=="WHITE" & cove$ETHNIC!="HISPANIC OR LATINO", 0, 1)
cove$aht<-as.numeric(difftime(strata_dt, dmy("23SEP2021"), units="days"))

# Obtain boost date, boost date + 13 days, and follow-up time
cove$Xstart_early<-as.numeric(difftime(ymd(cove$TR03SDT), dmy("23SEP2021"), units="days"))
cove$Xend_early<-as.numeric(difftime(ymd(cove$TR03SDT), dmy("23SEP2021"), units="days"))+13
cove$Xstart_full<-as.numeric(difftime(ymd(cove$TR03SDT), dmy("23SEP2021"), units="days"))+13
cove$futime<-as.numeric(difftime(cove$END_FU, dmy("23SEP2021"), units="days"))

# Obtain time to Delta/Omicron infection
cove <- cove %>% 
  mutate(
    Omicron_time = SARStte,
    Delta_time = SARStte,
    Omicron_status = ifelse(SARSevent==1 & LGROUP=="Omicron", 1, 0),
    Delta_status = ifelse(SARSevent==1 & LGROUP=="Delta", 1, 0),
  )

# Initialize `start`/`stop` data 
l0<-tmerge(data1=cove %>% 
             dplyr::select(SUBJID, STRATAR, SEX, MinorityInd, standardized_risk_score) %>% 
             distinct(), 
           data2=cove, id=SUBJID, tstop=futime)

# Create event for Omicron and Delta infections
l0<-tmerge(l0, cove, id=SUBJID, Omicron = event(Omicron_time, Omicron_status))
l0<-tmerge(l0, cove, id=SUBJID, Delta = event(Delta_time, Delta_status))

#Create time-dependent covariates for immediate vaccination and vaccination + 13 days
l0<-tmerge(l0, cove, id=SUBJID, early_boost_1 = tdc(Xstart_early))
l0<-tmerge(l0, cove, id=SUBJID, early_boost_2 = tdc(Xend_early))

# create variables for 'early boost' (within 13 days of vaccination) and 'late boost' (≥14 days of vaccination)
l0<-l0 %>% 
  mutate(early_boost = early_boost_1 - early_boost_2) 
l0<-l0 %>% 
  mutate(late_boost = ifelse(early_boost_1==1 & early_boost_2==1, 1, 0)) %>% 
  dplyr::select(-c("early_boost_1", "early_boost_2"))

l0<-left_join(l0, cove %>% dplyr::select(SUBJID, Xstart_early, Xstart_full, SARStte), by="SUBJID")

# Filter for observations before SARS-CoV-2 event
cove<-l0 %>% filter(tstop <= SARStte) 

### Delta endpoint

cove$vacc_time<-ifelse(is.na(cove$Xstart_early), Inf, cove$Xstart_early)
cove$vacc_status <- pmax(cove$early_boost, cove$late_boost)

omicron_events <- l0 %>% filter(Omicron==1) %>% mutate(tau = pmax(0, SARStte - Xstart_full, na.rm=TRUE), type="Omicron")
delta_events <- l0 %>% filter(Delta==1) %>% mutate(tau = pmax(0, SARStte - Xstart_full, na.rm=TRUE), type="Delta")

########################
# Create dataset for biofire data
########################

### Process biofire data

ag<-biofire %>% 
  dplyr::select(SUBJID, MBTEST, MBSCAT, mbdt_positive)

MBTESTs = biofire$MBTEST %>% unique()

id <- oracle %>% filter(PARAMCD=="TTCVDC1") %>% .$SUBJID %>% unique()

grid = expand.grid("SUBJID"=id, "MBTEST"=MBTESTs, "MBSCAT"="Nasopharyngeal Swab", "mbdt_positive"="")

ag<-rbind(ag,grid)

ag<-distinct(ag)

#2) join Part A, B, C Start and End dates to PE data

ag <- left_join(
  oracle %>% dplyr::select(SUBJID, START_A, END_A, START_OL, START_BOOST, END_FU, AP02SDT, AP02EDT, AP03SDT, AP03EDT, TR03SDT, EOSDT, EFFCODT, TRT01P, PPROTFL, FASFL, MITTFL, STRATAR, SEX, RACE, ETHNIC, TRTSEQA, PPPSFL,standardized_risk_score, ends_with("_FU")), 
  ag, by="SUBJID")

# Collect indicator and time-to-off-target events

boost_a <- ag %>% 
  filter(SUBJID %in% analysis_id) %>%
  mutate(
    `BOOST_nonsars` = case_when(
      END_FU <= ref_dt ~ NA,
      END_FU > ref_dt & is.na(dmy(mbdt_positive)) ~ 0,
      END_FU > ref_dt & !is.na(dmy(mbdt_positive)) ~ as.numeric( dplyr::between(dmy(mbdt_positive), ref_dt, pmin(END_FU, cens_dt)))),
    `tte_BOOST_nonsars` = case_when(
      END_FU <= ref_dt ~ NA,
      END_FU > ref_dt & is.na(dmy(mbdt_positive)) ~ as.numeric(difftime(pmin(END_FU, cens_dt), ref_dt, units="days")),
      END_FU > ref_dt & !is.na(dmy(mbdt_positive)) & dmy(mbdt_positive) <= ref_dt ~ as.numeric(difftime(pmin(END_FU, cens_dt), ref_dt, units="days")),
      END_FU > ref_dt & !is.na(dmy(mbdt_positive)) & dmy(mbdt_positive) > ref_dt ~ as.numeric(difftime(pmin(END_FU, dmy(mbdt_positive), cens_dt), ref_dt, units="days"))
    )
  )

# Filter for distinct entries
boost_a <- boost_a %>%
  dplyr::select(SUBJID, SEX, RACE, ETHNIC, TRTSEQA, PPROTFL, PPPSFL, MITTFL, STRATAR, standardized_risk_score, TR03SDT, tte_BOOST_nonsars, BOOST_nonsars, END_FU) %>% 
  group_by(SUBJID) %>% 
  distinct() %>%
  ungroup() %>%
  mutate(entry=0,
         URM_status = ifelse(RACE=="WHITE" & ETHNIC!="HISPANIC OR LATINO", 0, 1))

# Flag early boost date, late boost date, and FUtime, and at home testing date
boost_a$Xstart_early<-as.numeric(difftime(dmy(boost_a$TR03SDT), dmy("23SEP2021"), units="days"))
boost_a$Xend_early<-as.numeric(difftime(dmy(boost_a$TR03SDT), dmy("23SEP2021"), units="days"))+13
boost_a$Xstart_full<-as.numeric(difftime(dmy(boost_a$TR03SDT), dmy("23SEP2021"), units="days"))+13
boost_a$futime<-as.numeric(difftime(boost_a$END_FU, dmy("23SEP2021"), units="days"))
ath_thresh <- as.numeric(difftime(dmy("15JAN2022"), ref_dt, units="days"))

# Create dataset of non-sars-cov-2 (off-target) event indicators and event times (multiple records per PPT)
NS_events<-boost_a %>% dplyr::select(SUBJID, tte_BOOST_nonsars, BOOST_nonsars)

# Obtain vaccination times for each ppt
vaxx_times <- boost_a %>% 
  dplyr::select(SUBJID, 
                Xstart_early, 
                Xend_early) %>% 
  distinct()

# Join events and vaccination times
l<-left_join(NS_events %>% group_by(SUBJID) %>% mutate(futime = max(tte_BOOST_nonsars)), 
             vaxx_times, 
             by="SUBJID") %>% 
  arrange(SUBJID, tte_BOOST_nonsars) %>% 
  group_by(SUBJID) %>% 
  mutate(event_num=row_number()) %>% 
  pivot_wider(names_from=event_num, values_from=c(tte_BOOST_nonsars, BOOST_nonsars))

l0<-tmerge(data1=boost_a %>% dplyr::select(SUBJID, STRATAR, SEX, URM_status, standardized_risk_score) %>% distinct(), data2=l, id=SUBJID, tstop=futime)

l0<-tmerge(l0, l, id=SUBJID, infect = event(tte_BOOST_nonsars_1, BOOST_nonsars_1))
l0<-tmerge(l0, l, id=SUBJID, infect = event(tte_BOOST_nonsars_2, BOOST_nonsars_2))
l0<-tmerge(l0, l, id=SUBJID, infect = event(tte_BOOST_nonsars_3, BOOST_nonsars_3))
l0<-tmerge(l0, l, id=SUBJID, infect = event(tte_BOOST_nonsars_4, BOOST_nonsars_4))
l0<-tmerge(l0, l, id=SUBJID, early_boost_1 = tdc(Xstart_early))
l0<-tmerge(l0, l, id=SUBJID, early_boost_2 = tdc(Xend_early))
l0<-l0 %>% mutate(early_boost = early_boost_1 - early_boost_2) #%>% dplyr::select(-c("early_boost_1", "early_boost_2"))
l0<-l0 %>% mutate(late_boost = ifelse(early_boost_1==1 & early_boost_2==1, 1, 0)) %>% dplyr::select(-c("early_boost_1", "early_boost_2"))
l0<-left_join(l0, vaxx_times, by="SUBJID")

#replace unboosted with Infinite boost date
l0$standardized_risk_score<-ifelse(is.na(l0$standardized_risk_score), 0, l0$standardized_risk_score)
#l0$Xend_early<-ifelse(is.na(l0$Xend_early), Inf, l0$Xend_early)
l0$vacc_time<-ifelse(is.na(l0$Xstart_early), Inf, l0$Xstart_early)

nonsars_events <- l0 %>% 
  filter(infect ==1) %>% 
  mutate(tau = pmax(0, tstop - Xend_early, na.rm=TRUE), type="Non-SARS")

events <- bind_rows(
  delta_events %>% dplyr::select(1:7, 10, 11, 12, 13, 16),
  omicron_events %>% dplyr::select(1:7, 10, 11, 12, 13, 16),
  nonsars_events %>% dplyr::select(1:7, 9, 10, 11, 12, 15) %>% rename(MinorityInd = URM_status, Xstart_full=Xend_early)
)

# Select earliest event

events <- events %>% 
  group_by(SUBJID) %>% 
  slice(which.min(tstop)) %>% 
  ungroup()


#############################################################################################
# Period-specific Estimates
#############################################################################################

delta_period <- as.numeric(c(difftime(dmy("23SEP2021"), dmy("23SEP2021"), units="days"), difftime(dmy("11DEC2021"), dmy("23SEP2021"))))
omicron_period <- as.numeric(c(difftime(dmy("18DEC2021"), dmy("23SEP2021")), difftime(dmy("05APR2022"), dmy("23SEP2021"))))

#####################
# Delta period
#####################

#########
# Standard Cox regression
#########

delta_cove <- cove %>% 
  mutate(
    event = case_when(
      tstop <= delta_period[2] ~ Omicron + Delta,
      tstop > delta_period[2] ~ 0
    ),
    tstart = pmin(delta_period[2], tstart),
    tstop = pmin(delta_period[2], tstop)
  ) %>%
  filter(tstart < tstop)
  
fit_delta_period <- coxph(formula = Surv(time=tstart, time2=tstop, event=event) ~ 
        early_boost + 
        late_boost + 
        tt(vacc_time) + 
        strata(STRATAR) + 
        strata(SEX) + 
        strata(MinorityInd) + 
        standardized_risk_score, 
      tt = function(vacc_time, t, ...){
        # Only contributes if at least 14 days after vaccine
        as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
      }, 
      data=delta_cove)

max_tau_delta <- max(delta_cove$tstop[delta_cove$event==1] - delta_cove$Xstart_early[delta_cove$event==1], na.rm=TRUE)

grid <- cbind(seq(0, max_tau_delta)<14, seq(0, max_tau_delta)>=14, as.numeric(seq(0, max_tau_delta)>=14) * seq(0, max_tau_delta))

f <- grid %*% coef(fit_delta_period)[1:3]

ses <- sqrt(diag(grid %*% vcov(fit_delta_period)[1:3, 1:3] %*% t(grid)))

f_lci <- f - qnorm(0.975) * ses
f_uci <- f + qnorm(0.975) * ses

delta_out_cox<-data.frame(
  t = seq(0, max_tau_delta, by=1),
  type = rep("Cox PH", each=length(seq(0, max_tau_delta, by=1))),
  est= pmax(1-exp(f), 0),
  LCI=pmax(1-exp(f_uci), 0), 
  UCI=pmin(1-exp(f_lci), 1))

#########
# Our SLR approach
########

# events <- bind_rows(
#   delta_events %>% dplyr::select(1:7, 10, 11, 12, 13, 16),
#   omicron_events %>% dplyr::select(1:7, 10, 11, 12, 13, 16),
#   nonsars_events %>% dplyr::select(1:7, 10, 11, 12, 13, 16) %>% rename(MinorityInd = URM_status, Xstart_full=Xend_early)
# )

comp_risk=events %>% 
  dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
  mutate(J = case_when(
    type=="Omicron" ~ 2,
    type=="Delta" ~ 1,
    type=="Non-SARS" ~ 0
  )) %>%
  rename(`T`=tstop, 
         V=Xstart_full, 
         V_early = Xstart_early) %>%
  mutate(
    V = ifelse(is.na(V), Inf, V),
    V_early = ifelse(is.na(V_early), Inf, V_early)
  )

comp_risk_delta_period <- comp_risk %>%
  filter(`T` >= delta_period[1] & `T` <= delta_period[2]) %>%
  mutate(
    J = ifelse(J==0, 0, 1)
  )

# Sieve and TMLE estimators

sieve_delta <- sieve_partially_linear_logistic(comp_risk_delta_period,
                                               psi_delta = psi_d2_early,
                                               V.early.name = "V_early",
                                               monotone = FALSE,
                                               verbose = TRUE
)

tmle_delta <- tmle_iterative(comp_risk_delta_period,
                             psi_delta = psi_d2_early,
                             V.early.name = "V_early",
                             monotone = FALSE,
                             verbose = TRUE,
                             tol=1E-8)


# sieve_delta <- sieve_partially_linear_logistic(comp_risk_delta_period, 
#                                                psi_delta = psi_bs_early, 
#                                                first_nonboundary=0.6,
#                                                V.early.name = "V_early", 
#                                                df=6, 
#                                                monotone = FALSE,
#                                                intercept=TRUE,
#                                                verbose = TRUE
#                                                )
# 
# tmle_delta <- tmle_iterative(comp_risk_delta_period, 
#                              psi_delta = psi_bs_early, 
#                              first_nonboundary=0.6,
#                              V.early.name = "V_early", 
#                              df=6, 
#                              intercept=TRUE,
#                              monotone = FALSE,
#                              verbose = TRUE,
#                              tol=1E-8)

# Apply monotone corrections

#max tau is wrt 14 days post-boost
max_tau_delta <- max(comp_risk_delta_period$T - comp_risk_delta_period$V)

grid<-cbind(
  c(rep(1, 14), rep(0, max_tau_delta+1)),
  #c(rep(0, 14), rep(1, max_tau_delta+1)),
  rbind(matrix(0, nrow=14, ncol=ifelse(is.null(ncol(sieve_delta$basis)), 1, ncol(sieve_delta$basis))),
        predict(sieve_delta$basis, seq(0, max_tau_delta, by=1)))
)

grid<-cbind(c(rep(1, sum(seq(0, max_tau_delta+14, by=1) < 14)), rep(0, sum(seq(0, max_tau_delta+14, by=1)>=14))),
            c(rep(0, sum(seq(0, max_tau_delta+14, by=1)< 14)), rep(1, sum(seq(0, max_tau_delta+14, by=1)>=14))),
            c(rep(0, sum(seq(0, max_tau_delta+14, by=1)< 14)), rep(1, sum(seq(0, max_tau_delta+14, by=1)>=14))) *
              predict(sieve_delta$basis, seq(0, max_tau_delta+14, by=1)))

Amat = matrix(0, nrow = length(which(grid[,1]==0))-1, ncol= length(which(grid[,1]==0)))

for(i in 1:nrow(Amat)){
  
  Amat[i,i] <- -1
  Amat[i,(i+1)] <- 1
  
}

f_delta <- isotone_f(sieve_delta$beta_unconstr, 
          vcov = sieve_delta$cov, 
          grid = grid, 
          indices_to_monotonize = which(grid[,1]==0),
          Amat = Amat
          )

f_tmle <- isotone_f(tmle_delta$beta_unconstr, 
                     vcov = tmle_delta$cov, 
                     grid = grid, 
                     indices_to_monotonize = which(grid[,1]==0),
                     Amat = Amat
)

delta_out<-data.frame(
  t = rep(seq(0, max_tau_delta+14, by=1), 2),
  type = rep(c("Sieve", "TMLE"), each=length(seq(0, max_tau_delta+14, by=1))),
  est= pmax(1-exp(c(f_delta$f_mono, f_tmle$f_mono)), 0),
  LCI=pmax(1-exp(c(f_delta$f_uci, f_tmle$f_uci)), 0), 
  UCI=pmin(1-exp(c(f_delta$f_lci, f_tmle$f_lci)), 1))

p1<-bind_rows(
  delta_out_cox,
  delta_out
) %>%
  mutate(LCI = pmax(0, LCI)) %>%
  ggplot()+
  geom_jitter(aes(x=tau, y=0, color=factor(J)), 
              data=comp_risk_delta_period %>% 
                mutate(tau = pmax(0, T-V_early),
                       J = factor(J, levels=c(0,1), labels=c("Vaccine-Irrelevant ARI", "COVID-19"))), 
              alpha=0.25, height = 0.1, width=0)+
  scale_color_manual(values=c("black", "red"))+
  guides(color = guide_legend(title = NULL))+
  new_scale_color()+
  geom_line(aes(x=t, y=est, color=type), linewidth=1.35)+
  geom_ribbon(aes(x=t, ymin=pmax(0,LCI), ymax=UCI, fill=type), alpha=0.3)+
  scale_color_manual(values = c("grey", "orange", "#DAB1DA"))+
  scale_fill_manual(values = c("grey", "orange", "#DAB1DA"))+
  theme_bw()+
  labs(y="VE (1-HR)", x="Days since vaccination", title="Delta period (Sep 23, 2021-Dec 11, 2021)")+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title=NULL))+
  #facet_wrap(~variant, scales="free")+
  scale_y_continuous(labels=scales::percent, limits=c(-0.1,1), breaks=seq(0, 1, 0.1))+
  theme(legend.position="bottom", 
        strip.text=element_text(size=20), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        legend.text = element_text(size=16),
        strip.background = element_blank())


#########
# Standard Cox regression
#########

omicron_cove <- cove %>% 
  mutate(
    tstart = pmin(pmax(tstart, omicron_period[1]), omicron_period[2]),
    event = case_when(
      tstop <= omicron_period[2] ~ Omicron + Delta,
      tstop > omicron_period[2] ~ 0
    ),
    tstop = pmin(omicron_period[2], tstop)
  ) %>%
  filter(tstart < tstop)

fit_omicron_period <- coxph(formula = Surv(time=tstart, time2=tstop, event=event) ~ 
                            early_boost + 
                            late_boost + 
                            tt(vacc_time) + 
                            strata(STRATAR) + 
                            strata(SEX) + 
                            strata(MinorityInd) + 
                            standardized_risk_score, 
                          tt = function(vacc_time, t, ...){
                            # Only contributes if at least 14 days after vaccine
                            as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
                          }, 
                          data=omicron_cove)

max_tau_omicron <- max(omicron_cove$tstop[omicron_cove$event==1] - omicron_cove$vacc_time[omicron_cove$event==1], na.rm=TRUE)

grid <- cbind(seq(0, max_tau_omicron)<14, seq(0, max_tau_omicron)>=14, as.numeric(seq(0, max_tau_omicron)>=14) * seq(0, max_tau_omicron))

f <- grid %*% coef(fit_omicron_period)[1:3]

ses <- sqrt(diag(grid %*% vcov(fit_omicron_period)[1:3, 1:3] %*% t(grid)))

f_lci <- f - qnorm(0.975) * ses
f_uci <- f + qnorm(0.975) * ses

omicron_out_cox<-data.frame(
  t = seq(0, max_tau_omicron, by=1),
  type = rep("Cox PH", each=length(seq(0, max_tau_omicron, by=1))),
  est= pmax(1-exp(f), 0),
  LCI=pmax(1-exp(f_uci), 0), 
  UCI=pmin(1-exp(f_lci), 1))

#########
# Our SLR approach
########

comp_risk=events %>% 
  dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
  mutate(J = case_when(
    type=="Omicron" ~ 2,
    type=="Delta" ~ 1,
    type=="Non-SARS" ~ 0
  )) %>%
  rename(`T`=tstop, 
         V=Xstart_full, 
         V_early = Xstart_early) %>%
  mutate(
    V = ifelse(is.na(V), Inf, V),
    V_early = ifelse(is.na(V_early), Inf, V_early)
  )

comp_risk_omicron_period <- comp_risk %>%
  filter(`T` >= omicron_period[1] & `T` <= omicron_period[2]) %>%
  mutate(
    J = ifelse(J==0, 0, 1)
  )

# Sieve and TMLE estimators

sieve_omicron <- sieve_partially_linear_logistic(comp_risk_omicron_period, 
                                                 psi_delta = psi_d2_early, 
                                                 V.early.name = "V_early", 
                                                 verbose = TRUE
)

tmle_omicron <- tmle_iterative(comp_risk_omicron_period, 
                               psi_delta = psi_d2_early, 
                               V.early.name = "V_early", 
                               verbose = TRUE,
                               tol=1E-7)


# sieve_omicron <- sieve_partially_linear_logistic(comp_risk_omicron_period, 
#                                                psi_delta = psi_bs_early, 
#                                                first_nonboundary=0.05,
#                                                V.early.name = "V_early", 
#                                                df=6, 
#                                                monotone = FALSE,
#                                                intercept=TRUE,
#                                                verbose = TRUE
# )
# 
# tmle_omicron <- tmle_iterative(comp_risk_omicron_period, 
#                              psi_delta = psi_bs_early, 
#                              first_nonboundary=0.05,
#                              V.early.name = "V_early", 
#                              df=6, 
#                              intercept=TRUE,
#                              monotone = FALSE,
#                              verbose = TRUE,
#                              tol=1E-7)

# Apply monotone corrections

#max tau is wrt 14 days post-boost
max_tau_omicron <- max(comp_risk_omicron_period$T - comp_risk_omicron_period$V)

step_size=0.25
times_pre = seq(0, 14-step_size, by=step_size)
times_post = seq(0, max_tau_omicron, by=step_size)
times = c(times_pre, times_post+14)

grid<-cbind(c(rep(1, length(times_pre)), rep(0, length(times_post))),
            c(rep(0, length(times_pre)), rep(1, length(times_post))),
            rbind(matrix(0, nrow=length(times_pre), ncol=ifelse(is.null(ncol(sieve_omicron$basis)), 1, ncol(sieve_omicron$basis))),
                  predict(sieve_omicron$basis, times_post)))

Amat = matrix(0, nrow = length(which(grid[,1]==0))-1, ncol= length(which(grid[,1]==0)))

for(i in 1:nrow(Amat)){
  
  Amat[i,i] <- -1
  Amat[i,(i+1)] <- 1
  
}

f_omicron <- isotone_f(sieve_omicron$beta_unconstr, 
                     vcov = sieve_omicron$cov, 
                     grid = grid, 
                     indices_to_monotonize = which(grid[,1]==0),
                     Amat = Amat
)

f_tmle <- isotone_f(tmle_omicron$beta_unconstr, 
                    vcov = tmle_omicron$cov, 
                    grid = grid, 
                    indices_to_monotonize = which(grid[,1]==0),
                    Amat = Amat
)

omicron_out<-data.frame(
  t = rep(times, 2),
  type = rep(c("Sieve", "TMLE"), each=length(times)),
  est= pmax(1-exp(c(f_omicron$f_mono, f_tmle$f_mono)), 0),
  LCI=pmax(1-exp(c(f_omicron$f_uci, f_tmle$f_uci)), 0), 
  UCI=pmin(1-exp(c(f_omicron$f_lci, f_tmle$f_lci)), 1))

p2 <- bind_rows(
  omicron_out_cox,
  omicron_out
) %>%
  mutate(LCI = pmax(0, LCI)) %>%
  ggplot()+
  geom_jitter(aes(x=tau, y=0, color=factor(J)), 
              data=comp_risk_omicron_period %>% 
                mutate(tau = pmax(0, T-V_early),
                       J = factor(J, levels=c(0,1), labels=c("Vaccine-Irrelevant ARI", "COVID-19"))), 
              alpha=0.25, height = 0.1, width=0)+
  scale_color_manual(values=c("black", "red"))+
  guides(color = guide_legend(title = NULL))+
  new_scale_color()+
  geom_line(aes(x=t, y=est, color=type), linewidth=1.35)+
  geom_ribbon(aes(x=t, ymin=pmax(0,LCI), ymax=UCI, fill=type), alpha=0.3)+
  scale_color_manual(values = c("grey", "orange", "#DAB1DA"))+
  scale_fill_manual(values = c("grey", "orange", "#DAB1DA"))+
  theme_bw()+
  labs(y="VE (1-HR)", x="Days since vaccination", title="Omicron period (Dec 18, 2021-Apr 5, 2022)")+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title=NULL))+
  #facet_wrap(~variant, scales="free")+
  scale_y_continuous(labels=scales::percent, limits=c(-0.1,1), breaks=seq(0, 1, 0.1))+
  theme(legend.position="bottom", 
        strip.text=element_text(size=20), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        legend.text = element_text(size=16),
        strip.background = element_blank())


pdf("~/Desktop/VE_fig_1_LIN.pdf", width=9, height=4.5)
p1+p2+ plot_layout(guides="collect") & theme(legend.position="bottom")
dev.off()

#########################################################################
# Competing risk analysis
#########################################################################

##
# Standard Delta Cox regression
## 

fit_DELTA<-coxph(formula = Surv(time=tstart, time2=tstop, event=Delta) ~ 
                   early_boost + 
                   late_boost + 
                   tt(vacc_time) + 
                   strata(STRATAR) + 
                   strata(SEX) + 
                   strata(MinorityInd) + 
                   standardized_risk_score, 
                 tt = function(vacc_time, t, ...){
                   # Only contributes if at least 14 days after vaccine
                   as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
                 }, 
                 data=cove)

max_tau_delta <- cove %>% 
  filter(Delta==1) %>% 
  summarise(tau = pmax(0, tstop-Xstart_full, na.rm=TRUE)) %>% max(.$tau)

step_size=0.25
times_pre = seq(0, 14-step_size, by=step_size)
times_post = seq(0, max_tau_delta, by=step_size)
times = c(times_pre, times_post+14)

grid = cbind(
  c(rep(1, length(times_pre)), rep(0, length(times_post))),
  c(rep(0, length(times_pre)), rep(1, length(times_post))),
  I(times >= 14) * (times-14))


p<-length(coef(fit_DELTA))-1
res_delta_cox<-data.frame(t=times)
res_delta_cox$est<-apply(grid, 1, function(x) (x) %*% coef(fit_DELTA)[1:p])
res_delta_cox$se <- c(sqrt(apply(grid, 1, function(x) (x) %*% vcov(fit_DELTA)[1:p,1:p] %*% x)))
res_delta_cox$type="Delta"

res_cox<-rbind(res_delta_cox)

##
# Standard Omicron Cox regression
## 

fit_OMICRON<-coxph(formula = Surv(time=tstart, time2=tstop, event=Omicron) ~ 
                     early_boost + 
                     late_boost + 
                     tt(vacc_time) + 
                     strata(STRATAR) + 
                     strata(SEX) + 
                     strata(MinorityInd) + 
                     standardized_risk_score, 
                   tt = function(vacc_time, t, ...){
                     # Only contributes if at least 14 days after vaccine
                     as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
                   }, 
                   data=cove)


max_tau_omicron <- cove %>% 
  filter(Omicron==1) %>% 
  summarise(tau = pmax(0, tstop-Xstart_full, na.rm=TRUE)) %>% max(.$tau)

step_size=0.25
times_pre = seq(0, 14-step_size, by=step_size)
times_post = seq(0, max_tau_omicron, by=step_size)
times = c(times_pre, times_post+14)

grid = cbind(
  c(rep(1, length(times_pre)), rep(0, length(times_post))),
  c(rep(0, length(times_pre)), rep(1, length(times_post))),
  I(times >= 14) * (times-14))


p<-length(coef(fit_OMICRON))-1
res_omicron_cox<-data.frame(t=times)
res_omicron_cox$est<-apply(grid, 1, function(x) (x) %*% coef(fit_OMICRON)[1:p])
res_omicron_cox$se <- c(sqrt(apply(grid, 1, function(x) (x) %*% vcov(fit_OMICRON)[1:p,1:p] %*% x)))
res_omicron_cox$type="Omicron"

res_cox<-rbind(res_cox, res_omicron_cox)

##
# SLR
##

# GISAID data plot

GISAID_boost %>% 
  ggplot(aes(x=Week.prior.to, y=X..per.Country.and.Week, color=Value, group=Value))+
  geom_line()+
  theme(legend.position="none")

to_profile <- GISAID_boost %>% 
  filter(Value %in% c(
    "Former VOC Delta GK (B.1.617.2+AY.*) first detected in India",
    "Former VOC Omicron GRA (B.1.1.529+BA.*) first detected in Botswana/Hong Kong/South Africa")) %>%
  group_by(Week.prior.to) %>%
  mutate(
    Total_count = sum(Submission.Count),
  ) %>%
  group_by(Week.prior.to, Value) %>%
  summarise(
    prop = Submission.Count/Total_count
  ) %>%
  mutate(
    start = as.numeric(difftime(as.Date(Week.prior.to) - days(7), dmy("23SEP2021"))),
    end = as.numeric(difftime(as.Date(Week.prior.to) - days(1), dmy("23SEP2021")))
  )

to_profile <- to_profile %>% 
  mutate(Value = ifelse(grepl("Delta", Value), "Delta", "Omicron")) %>%
  pivot_wider(names_from=Value, values_from=prop) %>%
  mutate(Delta = ifelse(is.na(Delta), 0, Delta),
         Omicron = ifelse(is.na(Omicron), 0, Omicron))

fit_GISAID<-gam(Omicron ~ start, data=to_profile, family = binomial)

predict_GISAID<-function(t){
  
  preds<-predict(fit_GISAID, data.frame(start = t), type="response")
  cbind(1-preds,
        preds)
  
}

comp_risk=events %>% dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
  mutate(J = case_when(
    type=="Omicron" ~ 2,
    type=="Delta" ~ 1,
    type=="Non-SARS" ~ 0
  )) %>%
  rename(`T`=tstop, V=Xstart_full, V_early = Xstart_early) %>%
  mutate(
    V = ifelse(is.na(V), Inf, V),
    V_early = ifelse(is.na(V_early), Inf, V_early)
  )

comp_risk$p_Delta <- 1-predict(fit_GISAID, data.frame(start = comp_risk$T), type="response")
comp_risk$p_Omicron <- predict(fit_GISAID, data.frame(start = comp_risk$T), type="response")

source("~/Desktop/GitHub/TimeVaryingVE/sieve_multinomial_1.r")
source("~/Desktop/GitHub/TimeVaryingVE/multinomial_TMLE.r")

res_comprisk_sieve<-sieve_multinomial(dat=comp_risk, 
                                   V.early.name = "V_early", 
                                   psi_delta = psi_d2_early, 
                                   monotone = FALSE,
                                   #df=6,
                                   maxiter=50, 
                                   tol = 1E-6,
                                   calc_mixture = FALSE, 
                                   pk_names = "p_Omicron")

res_comprisk_TMLE <- tmle_multinomial(dat = comp_risk %>% mutate(
                                          S = J,
                                          J = ifelse(J>0, 1, 0)), 
                                      psi_delta = psi_d2_early, 
                                      m = 2, 
                                      p_s_fun = predict_GISAID, 
                                      kernel_shape = "epanechnikov", 
                                      automatic_bw = FALSE,
                                      verbose=TRUE)

# var_f = diag((grid[,2:4] %*% list_out$var_beta1) %*% t(grid[,2:4]))
# 
# out<-data.frame(t=grid[,1], 
#                 VE=grid[,2:4] %*% beta[1:3], 
#                 VE_lci = (grid[,2:4] %*% beta[1:3]) + qnorm(0.975) * sqrt(var_f),
#                 VE_uci = (grid[,2:4] %*% beta[1:3]) - qnorm(0.975) * sqrt(var_f))
# 
# out$VE<-1-exp(out$VE)
# out$VE_lci<-1-exp(out$VE_lci)
# out$VE_uci<-1-exp(out$VE_uci)

# Collect results

max_tau_delta <- max(pmax(comp_risk$T[comp_risk$J==1] - comp_risk$V[comp_risk$J==1], 0), na.rm=TRUE)
max_tau_omicron <- max(pmax(comp_risk$T[comp_risk$J==2] - comp_risk$V[comp_risk$J==2], 0), na.rm=TRUE)

step_size=0.25
times_pre = seq(0, 14-step_size, by=step_size)
times_post = seq(0, max_tau_delta, by=step_size)
times = c(times_pre, times_post+14)

grid_delta<-cbind(
  c(rep(1, length(times_pre)), rep(0, length(times_post))),
  c(rep(0, length(times_pre)), rep(1, length(times_post))),
  rbind(matrix(0, nrow=length(times_pre), ncol=1),
        predict(res_comprisk_sieve$basis, times_post))
)

Amat = matrix(0, nrow = length(which(grid_delta[,1]==0))-1, ncol= length(which(grid_delta[,1]==0)))

for(i in 1:nrow(Amat)){
  
  Amat[i,i] <- -1
  Amat[i,(i+1)] <- 1
  
}

f_delta_sieve <- isotone_f(res_comprisk_sieve$beta_1, 
                       vcov = res_comprisk_sieve$cov[1:3, 1:3], 
                       grid = grid_delta, 
                       indices_to_monotonize = which(grid_delta[,1]==0),
                       Amat = Amat
)

f_delta_tmle <- isotone_f(res_comprisk_TMLE$beta[1:3], 
                     vcov = res_comprisk_TMLE$var_beta1, 
                     grid = grid_delta, 
                     indices_to_monotonize = which(grid_delta[,1]==0),
                     Amat = Amat
)

delta_out<-data.frame(
  t = rep(times, 2),
  type = rep(c("Sieve", "TMLE"), each=length(times)),
  est= pmax(1-exp(c(f_delta_sieve$f_mono, f_delta_tmle$f_mono)), 0),
  LCI=pmax(1-exp(c(f_delta_sieve$f_uci, f_delta_tmle$f_uci)), 0), 
  UCI=pmin(1-exp(c(f_delta_sieve$f_lci, f_delta_tmle$f_lci)), 1))

step_size=0.25
times_pre = seq(0, 14-step_size, by=step_size)
times_post = seq(0, max_tau_omicron, by=step_size)
times = c(times_pre, times_post+14)

grid_omicron<-cbind(
  c(rep(1, length(times_pre)), rep(0, length(times_post))),
  c(rep(0, length(times_pre)), rep(1, length(times_post))),
  rbind(matrix(0, nrow=length(times_pre), ncol=1),
        predict(res_comprisk_sieve$basis, times_post))
)

Amat = matrix(0, nrow = length(which(grid_omicron[,1]==0))-1, ncol= length(which(grid_omicron[,1]==0)))

for(i in 1:nrow(Amat)){
  
  Amat[i,i] <- -1
  Amat[i,(i+1)] <- 1
  
}

f_omicron_sieve <- isotone_f(res_comprisk_sieve$beta_2, 
                           vcov = res_comprisk_sieve$cov[4:6, 4:6], 
                           grid = grid_omicron, 
                           indices_to_monotonize = which(grid_omicron[,1]==0),
                           Amat = Amat
)

f_omicron_tmle <- isotone_f(res_comprisk_TMLE$beta[4:6], 
                          vcov = res_comprisk_TMLE$var_beta2, 
                          grid = grid_omicron, 
                          indices_to_monotonize = which(grid_omicron[,1]==0),
                          Amat = Amat
)

omicron_out<-data.frame(
  t = rep(times, 2),
  type = rep(c("Sieve", "TMLE"), each=length(times)),
  est= 1-exp(c(f_omicron_sieve$f_mono, f_omicron_tmle$f_mono)),
  LCI= 1-exp(c(f_omicron_sieve$f_uci, f_omicron_tmle$f_uci)), 
  UCI= 1-exp(c(f_omicron_sieve$f_lci, f_omicron_tmle$f_lci)))

p1 <- bind_rows(
  res_delta_cox %>% 
    mutate(type="Cox PH") %>% 
    mutate(f_uci = est + qnorm(0.975) * se,
           f_lci = est - qnorm(0.975) * se) %>% 
    mutate(est = 1-exp(est), LCI = 1-exp(f_uci), UCI = 1-exp(f_lci)),
  delta_out
) %>% 
  mutate(LCI = pmax(0, LCI)) %>%
  ggplot()+
  geom_jitter(aes(x=tau, y=0, color=factor(J)), 
              data=comp_risk %>% 
                filter(J!=2) %>%
                mutate(tau = pmax(0, T-V_early),
                       max_tau_delta = max(pmax(0, T-V_early) * I(J==1)),
                       J = factor(J, levels=c(0,1), labels=c("Vaccine-Irrelevant ARI", "COVID-19"))) %>%
                filter(tau <= max_tau_delta), 
              alpha=0.25, height = 0.1, width=0)+
  scale_color_manual(values=c("black", "red"))+
  guides(color = guide_legend(title = NULL, order=1))+
  new_scale_color()+
  geom_line(aes(x=t, y=est, color=type), linewidth=1.35)+
  geom_ribbon(aes(x=t, ymin=pmax(0,LCI), ymax=UCI, fill=type), alpha=0.3)+
  scale_color_manual(values = c("grey", "orange", "#DAB1DA"))+
  scale_fill_manual(values = c("grey", "orange", "#DAB1DA"))+
  theme_bw()+
  labs(y="VE (1-HR)", x="Days since vaccination", title="Delta COVID-19")+
  guides(color = guide_legend(title = NULL, order=2), fill = guide_legend(title=NULL, order=2))+
  #facet_wrap(~variant, scales="free")+
  scale_y_continuous(labels=scales::percent, limits=c(-0.1,1), breaks=seq(0, 1, 0.1))+
  theme(legend.position="bottom", 
        strip.text=element_text(size=20), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        legend.text = element_text(size=16),
        strip.background = element_blank())


p2 <- bind_rows(
  res_omicron_cox %>% 
    mutate(type="Cox PH") %>% 
    mutate(f_uci = est + qnorm(0.975) * se,
           f_lci = est - qnorm(0.975) * se) %>% 
    mutate(est = 1-exp(est), LCI = 1-exp(f_uci), UCI = 1-exp(f_lci)),
  omicron_out
) %>% 
  mutate(LCI = pmax(0, LCI)) %>%
  ggplot()+
  geom_jitter(aes(x=tau, y=0, color=factor(J)), 
              data=comp_risk %>% 
                filter(J!=1) %>%
                mutate(tau = pmax(0, T-V_early),
                       max_tau_omicron = max(pmax(0, T-V_early) * I(J==2)),
                       J = factor(J, levels=c(0,2), labels=c("Vaccine-Irrelevant ARI", "COVID-19"))) %>%
                filter(tau <= max_tau_omicron), 
              alpha=0.25, height = 0.1, width=0)+
  scale_color_manual(values=c("black", "red"))+
  guides(color = guide_legend(title = NULL, order=1))+
  new_scale_color()+
  geom_line(aes(x=t, y=est, color=type), linewidth=1.35)+
  geom_ribbon(aes(x=t, ymin=pmax(0,LCI), ymax=UCI, fill=type), alpha=0.3)+
  scale_color_manual(values = c("grey", "orange", "#DAB1DA"))+
  scale_fill_manual(values = c("grey", "orange", "#DAB1DA"))+
  theme_bw()+
  labs(y="VE (1-HR)", x="Days since vaccination", title="Omicron COVID-19")+
  guides(color = guide_legend(title = NULL, order=2), fill = guide_legend(title=NULL, order=2))+
  #facet_wrap(~variant, scales="free")+
  scale_y_continuous(labels=scales::percent, limits=c(-0.1,1), breaks=seq(0, 1, 0.1))+
  theme(legend.position="bottom", 
        strip.text=element_text(size=20), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        legend.text = element_text(size=16),
        strip.background = element_blank())


pdf("~/Desktop/VE_fig_2_LIN.pdf", width=9, height=4.5)
p1+p2+ plot_layout(guides="collect") & theme(legend.position="bottom")
dev.off()

# ####################################################################################
# # Extra
# ####################################################################################
# 
# ########
# # Process results of model fit
# ########
# 
# max_tau_delta <- cove %>% filter(Delta==1) %>% summarise(tau = pmax(0, tstop-vacc_time, na.rm=TRUE)) %>% quantile(.$tau, 0.99)
# 
# event_time <- seq(0, max_tau_delta)
# 
# n_early <- sum(event_time<14)
# 
# times_lin = cbind(
#   c(rep(1, n_early), rep(0, length(event_time)-n_early)),
#   c(rep(0, n_early), rep(1, length(event_time)-n_early)), event_time)
# 
# res_delta_cox<-data.frame(times=times_lin[,3])
# 
# times_lin[,3]<-ifelse(times_lin[,3]<14, 0, times_lin[,3])
# 
# p<-length(coef(fit_DELTA))-1
# 
# res_delta_cox$est<-apply(times_lin, 1, function(x) (x) %*% coef(fit_DELTA)[1:p])
# res_delta_cox$se <- c(sqrt(apply(times_lin, 1, function(x) (x) %*% vcov(fit_DELTA)[1:p,1:p] %*% x)))
# res_delta_cox$type="Delta"
# 
# res_cox<-rbind(res_delta_cox)
# 
# ########################
# # Fit model to Omicron SARS-CoV-2 -- linear
# #########################
# 
# fit_OMICRON<-coxph(formula = Surv(time=tstart, time2=tstop, event=Omicron) ~ 
#                    early_boost + 
#                    late_boost + 
#                    tt(vacc_time) + 
#                    strata(STRATAR) + 
#                    strata(SEX) + 
#                    strata(MinorityInd) + 
#                    standardized_risk_score, 
#                  tt = function(vacc_time, t, ...){
#                    # Only contributes if at least 14 days after vaccine
#                    as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
#                  }, 
#                  data=cove)
# 
# ########
# # Process results of model fit
# ########
# 
# max_tau_omicron <- cove %>% filter(Delta==1) %>% summarise(tau = pmax(0, tstop-vacc_time, na.rm=TRUE)) %>% quantile(.$tau, 0.99)
# 
# event_time <- seq(0, max_tau_omicron)
# 
# n_early <- sum(event_time<14)
# 
# times_lin = cbind(
#   c(rep(1, n_early), rep(0, length(event_time)-n_early)),
#   c(rep(0, n_early), rep(1, length(event_time)-n_early)), event_time)
# 
# res_omicron_cox<-data.frame(times=times_lin[,3])
# 
# times_lin[,3]<-ifelse(times_lin[,3]<14, 0, times_lin[,3])
# 
# p<-length(coef(fit_OMICRON))-1
# 
# res_omicron_cox$est<-apply(times_lin, 1, function(x) (x) %*% coef(fit_OMICRON)[1:p])
# res_omicron_cox$se <- c(sqrt(apply(times_lin, 1, function(x) (x) %*% vcov(fit_OMICRON)[1:p,1:p] %*% x)))
# res_omicron_cox$type="Omicron"
# 
# res_cox<-rbind(res_cox, res_omicron_cox)
# 
# res_cox$VE = 1-exp(res_cox$est)
# res_cox$VE_lci = 1-exp(res_cox$est + qnorm(0.975) * res_cox$se)
# res_cox$VE_uci = 1-exp(res_cox$est - qnorm(0.975) * res_cox$se)
# 
# ########################
# # Fit model to Omicron SARS-CoV-2 -- bs
# #########################
# 
# # fit_OMICRON<-coxph(formula = Surv(time=tstart, time2=tstop, event=Omicron) ~ 
# #                      early_boost + 
# #                      #late_boost + 
# #                      tt(vacc_time) +
# #                      #strata(STRATAR) + 
# #                      #strata(SEX) + 
# #                      #strata(MinorityInd) + 
# #                      standardized_risk_score, 
# #                    tt = function(vacc_time, t, ...){
# #                      # Only contributes if at least 14 days after vaccine
# #                      #as.numeric(t-vacc_time >= 14) * pmax(14, t - vacc_time, na.rm=TRUE)
# #                      tau = pmax(14, t - vacc_time, na.rm=TRUE)
# #                      
# #                      K = 6 - 3
# #                      
# #                      cuts <- seq(1, K) / (K + 1)
# #                      
# #                      knots <- quantile(unique(tau), cuts)
# #                      
# #                      basis <- bs(tau, df=6, degree = 3, knots = knots, intercept=TRUE)
# #                      #basis <- apply(basis, MARGIN=2, FUN=function(x){x-mean(x)})
# #                      as.numeric(pmax(0, t - vacc_time, na.rm=TRUE) >= 14) * basis
# #                    },
# #                    data=cove)
# # 
# # tau = pmax(14, cove$tstop - cove$vacc_time, na.rm=TRUE)
# # 
# # K = 6 - 3
# # 
# # cuts <- seq(1, K) / (K + 1)
# # 
# # knots <- quantile(unique(tau), cuts)
# # 
# # basis <- bs(tau, df=6, degree = 3, knots = knots, intercept=TRUE)
# # 
# # max_tau_omicron = max(tau[cove$Omicron==1])
# # 
# # eval_times<-seq(0, max_tau_omicron, 1)
# # 
# # b_mat <- cbind(c(rep(1, n_early), rep(0, length(eval_times)-n_early)),
# #       #c(rep(0, n_early), rep(1, length(eval_times)-n_early)),
# #       as.numeric(eval_times >= 14) * predict(basis, eval_times))
# # 
# # # times_lin = cbind(
# # #   c(rep(1, n_early), rep(0, length(event_time)-n_early)),
# # #   c(rep(0, n_early), rep(1, length(event_time)-n_early)), event_time)
# # 
# # res_omicron_cox<-data.frame(times=eval_times)
# # 
# # p<-length(coef(fit_OMICRON))-1
# # 
# # res_omicron_cox$est<-apply(b_mat, 1, function(x) (x) %*% coef(fit_OMICRON)[1:p])
# # res_omicron_cox$se <- c(sqrt(apply(b_mat, 1, function(x) (x) %*% vcov(fit_OMICRON)[1:p,1:p] %*% x)))
# # res_omicron_cox$type="Omicron"
# # 
# # res_cox<-bind_rows(res_cox, res_omicron_cox)
# # 
# # res_cox$VE = 1-exp(res_cox$est)
# # res_cox$VE_lci = 1-exp(res_cox$est + qnorm(0.975) * res_cox$se)
# # res_cox$VE_uci = 1-exp(res_cox$est - qnorm(0.975) * res_cox$se)
# 
# ###################################################
# # Fit B-spline multinomial WITH GISAID profile out
# ###################################################
# 
# comp_risk=events %>% dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
#   mutate(J = case_when(
#     type=="Omicron" ~ 2,
#     type=="Delta" ~ 1,
#     type=="Non-SARS" ~ 0
#   )) %>%
#   rename(`T`=tstop, V=Xstart_full, V_early = Xstart_early) %>%
#   mutate(
#     V = ifelse(is.na(V), Inf, V),
#     V_early = ifelse(is.na(V_early), Inf, V_early)
#   )
# 
# # GISAID data plot
# GISAID_boost %>% 
#   ggplot(aes(x=Week.prior.to, y=X..per.Country.and.Week, color=Value, group=Value))+
#   geom_line()+
#   theme(legend.position="none")
# 
# to_profile <- GISAID_boost %>% 
#   filter(Value %in% c(
#     "Former VOC Delta GK (B.1.617.2+AY.*) first detected in India",
#     "Former VOC Omicron GRA (B.1.1.529+BA.*) first detected in Botswana/Hong Kong/South Africa")) %>%
#   group_by(Week.prior.to) %>%
#   mutate(
#     Total_count = sum(Submission.Count),
#   ) %>%
#   group_by(Week.prior.to, Value) %>%
#   summarise(
#     prop = Submission.Count/Total_count
#   ) %>%
#   mutate(
#     start = as.numeric(difftime(as.Date(Week.prior.to) - days(7), dmy("23SEP2021"))),
#     end = as.numeric(difftime(as.Date(Week.prior.to) - days(1), dmy("23SEP2021")))
#   )
# 
# to_profile <- to_profile %>% 
#   mutate(Value = ifelse(grepl("Delta", Value), "Delta", "Omicron")) %>%
#   pivot_wider(names_from=Value, values_from=prop) %>%
#   mutate(Delta = ifelse(is.na(Delta), 0, Delta),
#          Omicron = ifelse(is.na(Omicron), 0, Omicron))
# 
# fit<-gam(Omicron ~ start, data=to_profile, family = binomial)
# 
# # comp_risk$p_Delta <- NA
# # comp_risk$p_Omicron <- NA
# # for(j in 1:nrow(comp_risk)){
# #   
# #   indx<-which(between(comp_risk$`T`[j], to_profile$start, to_profile$end))
# #   #print(j)
# #   #print(indx)
# #   
# #   comp_risk$p_Delta[j]<-to_profile$Delta[indx]
# #   comp_risk$p_Omicron[j]<-to_profile$Omicron[indx]
# # }
# 
# comp_risk$p_Delta <- 1-predict(fit, data.frame(start = comp_risk$T), type="response")
# comp_risk$p_Omicron <- predict(fit, data.frame(start = comp_risk$T), type="response")
# 
# type="bs"
# 
# if(type=="bs"){
#   
#   Amat <- matrix(0, ncol=8, nrow=8)
#   
#   for(i in 1:nrow(Amat)){
#     
#     if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#       
#       Amat[i,i]<--1
#       
#     }else{ # Force VE to wane (cannot increase over time)
#       
#       Amat[i,(i-1)] <- -1
#       Amat[i,i] <- 1
#       
#     }
#     
#   }
#   
#   res_comprisk<-sieve_multinomial_EM(dat=comp_risk, 
#                                      V.early.name = "V_early", 
#                                      psi_delta = psi_bs_early, 
#                                      monotone = TRUE, 
#                                      Amat=Amat, 
#                                      df=6,
#                                      maxiter=50, 
#                                      tol = 1E-3,
#                                      calc_mixture = FALSE, 
#                                      pk_names = "p_Omicron")
#   
# }else{
#   Amat <- matrix(0, ncol=3, nrow=3)
#   
#   for(i in 1:nrow(Amat)){
#     
#     if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#       
#       Amat[i,i]<--1
#       
#     }else{ # Force VE to wane (cannot increase over time)
#       
#       Amat[i,i] <- 1
#       
#     }
#   }
#   
#   res_comprisk<-sieve_multinomial_EM(dat=comp_risk, 
#                                      V.early.name = "V_early", 
#                                      psi_delta = psi_d2_early, 
#                                      monotone = TRUE, 
#                                      Amat=Amat, 
#                                      maxiter=100, tol=1E-6)
#   
# }
# 
# grid <- cbind(as.numeric(seq(0, max_tau_delta) <= 13),
#               #as.numeric(seq(0, max_tau) >= 14),
#               as.numeric(seq(0, max_tau_delta) >= 14) * predict(res_comprisk$basis, seq(0, max_tau_delta)))
# 
# res_SLR_delta <- data.frame(t=seq(0, max_tau_delta))
# 
# res_SLR_delta$VE_sieve <- 1-exp(grid %*% res_comprisk$beta_1)
# 
# CI_sieve_delta <- monotone_CI_MC(beta = res_comprisk$beta_1, 
#                                  vcov = res_comprisk$cov[1:length(res_comprisk$beta_1), 1:length(res_comprisk$beta_1)], 
#                                  Amat = res_comprisk$Amat, 
#                                  w = 1/(res_comprisk$se[1:length(res_comprisk$beta_1)])^2,
#                                  grid = grid,
#                                  M=10000,
#                                  seed=47)
# 
# res_SLR_delta$VE_sieve_LCI <- 1-exp(CI_sieve_delta[,2])
# res_SLR_delta$VE_sieve_UCI <- 1-exp(CI_sieve_delta[,1])
# 
# grid <- cbind(as.numeric(seq(0, max_tau_omicron) <= 13),
#               #as.numeric(seq(0, max_tau) >= 14),
#               as.numeric(seq(0, max_tau_omicron) >= 14) * predict(res_comprisk$basis, seq(0, max_tau_omicron)))
# 
# res_SLR_omicron <- data.frame(t=seq(0, max_tau_omicron))
# 
# res_SLR_omicron$VE_sieve<- 1-exp(grid %*% res_comprisk$beta_2)
# 
# CI_sieve_omicron <- monotone_CI_MC(beta = res_comprisk$beta_2, 
#                                    vcov = res_comprisk$cov[(length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1)), 
#                                                            (length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1))], 
#                                    Amat = res_comprisk$Amat, 
#                                    w = 1/(res_comprisk$se[(length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1))])^2,
#                                    grid = grid,
#                                    M=10000,
#                                    seed=47)
# 
# res_SLR_omicron$VE_sieve_LCI <- 1-exp(CI_sieve_omicron[,2])
# res_SLR_omicron$VE_sieve_UCI <- 1-exp(CI_sieve_omicron[,1])
# 
# res_SLR <- bind_rows(
#   res_SLR_delta %>% mutate(variant="Delta"),
#   res_SLR_omicron %>% mutate(variant="Omicron")
# )
# 
# bind_rows(
#   res_SLR %>% 
#     rename(times = t, VE = VE_sieve, VE_lci = VE_sieve_LCI, VE_uci = VE_sieve_UCI) %>% mutate(method="SLR (Sieve)"),
#   res_cox %>%
#     dplyr::select(-c("est", "se")) %>% rename(variant = type) %>% mutate(method = "Cox PH") %>%
#     mutate(flag = case_when(
#       variant == "Delta" ~ times <= max_tau_delta,
#       variant == "Omicron" ~ times <= max_tau_omicron
#     )) %>%
#     filter(flag==TRUE) %>%
#     dplyr::select(-flag)
# ) %>%
#   mutate(VE_lci = pmax(0, VE_lci)) %>%
#   ggplot()+
#   geom_line(aes(x=times, y=VE, color=method), linewidth=1.35)+
#   geom_ribbon(aes(x=times, ymin=pmax(0,VE_lci), ymax=VE_uci, fill=method), alpha=0.3)+
#   scale_color_manual(values = c("grey", "orange"))+
#   scale_fill_manual(values = c("grey", "orange"))+
#   theme_bw()+
#   labs(y="VE (1-HR)", x="Days since vaccination", fill=NULL, color=NULL)+
#   facet_wrap(~variant, scales="free")+
#   scale_y_continuous(labels=scales::percent, limits=c(0,1))+
#   theme(legend.position="bottom", 
#         strip.text=element_text(size=20), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         legend.text = element_text(size=16),
#         strip.background = element_blank())
# 
# 
# ################################################################################################
# # FINAL B-spline multinomial
# ################################################################################################
# 
# comp_risk=events %>% dplyr::select(tstop, Xstart_early, Xstart_full, type) %>%
#   mutate(J = case_when(
#     type=="Omicron" ~ 2,
#     type=="Delta" ~ 1,
#     type=="Non-SARS" ~ 0
#   )) %>%
#   rename(`T`=tstop, V=Xstart_full, V_early = Xstart_early) %>%
#   mutate(
#     V = ifelse(is.na(V), Inf, V),
#     V_early = ifelse(is.na(V_early), Inf, V_early)
#   )
# 
# max_tau_delta <- comp_risk %>% filter(J==1) %>% summarise(tau = pmax(0, T-V)) %>% summarise(tau_delta = quantile(tau, 0.99)) %>% .$tau_delta
# max_tau_omicron <- comp_risk %>% filter(J==2) %>% summarise(tau = pmax(0, T-V)) %>% summarise(tau_omicron = quantile(tau, 0.99)) %>% .$tau_omicron
# 
# # DEPRECIATED
# # res_comprisk<-sieve_multinomial(dat=comp_risk, V.early.name = "V_early", psi_delta = psi_bs_early, df=6, monotone = TRUE, Amat=Amat)
# 
# type="bs"
# 
# if(type=="bs"){
# 
# Amat <- matrix(0, ncol=8, nrow=8)
# 
# for(i in 1:nrow(Amat)){
#   
#   if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#     
#     Amat[i,i]<--1
#     
#   }else{ # Force VE to wane (cannot increase over time)
#     
#     Amat[i,(i-1)] <- -1
#     Amat[i,i] <- 1
#     
#   }
#   
# }
# 
# res_comprisk<-sieve_multinomial_EM(dat=comp_risk, 
#                                    V.early.name = "V_early", 
#                                    psi_delta = psi_bs_early, 
#                                    monotone = TRUE, 
#                                    Amat=Amat, 
#                                    df=6,
#                                    maxiter=50, 
#                                    tol = 1E-3)
# 
# }else{
#   Amat <- matrix(0, ncol=3, nrow=3)
#   
#   for(i in 1:nrow(Amat)){
#     
#     if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#       
#       Amat[i,i]<--1
#       
#     }else{ # Force VE to wane (cannot increase over time)
#       
#       Amat[i,i] <- 1
#       
#     }
#   }
#   
#   res_comprisk<-sieve_multinomial_EM(dat=comp_risk, 
#                                      V.early.name = "V_early", 
#                                      psi_delta = psi_d2_early, 
#                                      monotone = TRUE, 
#                                      Amat=Amat, 
#                                      maxiter=100, tol=1E-6)
#   
# }
# 
# 
# grid <- cbind(as.numeric(seq(0, max_tau_delta) <= 13),
#               #as.numeric(seq(0, max_tau) >= 14),
#               as.numeric(seq(0, max_tau_delta) >= 14) * predict(res_comprisk$basis, seq(0, max_tau_delta)))
# 
# res_SLR_delta <- data.frame(t=seq(0, max_tau_delta))
# 
# res_SLR_delta$VE_sieve <- 1-exp(grid %*% res_comprisk$beta_1)
# 
# CI_sieve_delta <- monotone_CI_MC(beta = res_comprisk$beta_1, 
#                                  vcov = res_comprisk$cov[1:length(res_comprisk$beta_1), 1:length(res_comprisk$beta_1)], 
#                                  Amat = res_comprisk$Amat, 
#                                  w = 1/(res_comprisk$se[1:length(res_comprisk$beta_1)])^2,
#                                  grid = grid,
#                                  M=10000,
#                                  seed=47)
# 
# res_SLR_delta$VE_sieve_LCI <- 1-exp(CI_sieve_delta[,2])
# res_SLR_delta$VE_sieve_UCI <- 1-exp(CI_sieve_delta[,1])
# 
# grid <- cbind(as.numeric(seq(0, max_tau_omicron) <= 13),
#               #as.numeric(seq(0, max_tau) >= 14),
#               as.numeric(seq(0, max_tau_omicron) >= 14) * predict(res_comprisk$basis, seq(0, max_tau_omicron)))
# 
# res_SLR_omicron <- data.frame(t=seq(0, max_tau_omicron))
# 
# res_SLR_omicron$VE_sieve<- 1-exp(grid %*% res_comprisk$beta_2)
# 
# CI_sieve_omicron <- monotone_CI_MC(beta = res_comprisk$beta_2, 
#                                    vcov = res_comprisk$cov[(length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1)), 
#                                                            (length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1))], 
#                                    Amat = res_comprisk$Amat, 
#                                    w = 1/(res_comprisk$se[(length(res_comprisk$beta_1)+1):(2*length(res_comprisk$beta_1))])^2,
#                                    grid = grid,
#                                    M=10000,
#                                    seed=47)
# 
# res_SLR_omicron$VE_sieve_LCI <- 1-exp(CI_sieve_omicron[,2])
# res_SLR_omicron$VE_sieve_UCI <- 1-exp(CI_sieve_omicron[,1])
# 
# res_SLR <- bind_rows(
#   res_SLR_delta %>% mutate(variant="Delta"),
#   res_SLR_omicron %>% mutate(variant="Omicron")
# )
# 
# bind_rows(
#   res_SLR %>% 
#     rename(times = t, VE = VE_sieve, VE_lci = VE_sieve_LCI, VE_uci = VE_sieve_UCI) %>% mutate(method="SLR (Sieve)"),
#   res_cox %>%
#     dplyr::select(-c("est", "se")) %>% rename(variant = type) %>% mutate(method = "Cox PH") %>%
#     mutate(flag = case_when(
#       variant == "Delta" ~ times <= max_tau_delta,
#       variant == "Omicron" ~ times <= max_tau_omicron
#     )) %>%
#     filter(flag==TRUE) %>%
#     dplyr::select(-flag)
#   ) %>%
#   mutate(VE_lci = pmax(0, VE_lci)) %>%
#   ggplot()+
#   geom_line(aes(x=times, y=VE, color=method), linewidth=1.35)+
#   geom_ribbon(aes(x=times, ymin=pmax(0,VE_lci), ymax=VE_uci, fill=method), alpha=0.3)+
#   scale_color_manual(values = c("grey", "orange"))+
#   scale_fill_manual(values = c("grey", "orange"))+
#   theme_bw()+
#   labs(y="VE (1-HR)", x="Days since vaccination", fill=NULL, color=NULL)+
#   facet_wrap(~variant, scales="free")+
#   scale_y_continuous(labels=scales::percent, limits=c(0,1))+
#   theme(legend.position="bottom", 
#         strip.text=element_text(size=20), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         legend.text = element_text(size=16),
#         strip.background = element_blank())
# 
# 
# 
# ########################
# # Fit Sieve & TMLE models to Delta
# #########################
# 
# #max_t_delta <- as.numeric(difftime(dmy("15DEC2021"), dmy("23SEP2021")))
# max_t_delta <- cove %>% filter(Delta==1) %>% .$tstop %>% max()
# scale_f=365.25
# 
# delta_eligible <- events %>% 
#   filter(tstop <= max_t_delta & type != "Omicron") %>%
#   dplyr::select(Xstart_full, Xstart_early, tstop, type) %>%
#   mutate(tstop = tstop/scale_f, 
#          Xstart_early = Xstart_early/scale_f,
#          Xstart_full = Xstart_full/scale_f,
#          type = ifelse(type=="Delta", 1, 0)) %>%
#   mutate(
#     Xstart_early = case_when(
#       is.na(Xstart_early) ~ Inf,
#       .default = Xstart_early
#     ),
#     Xstart_full = case_when(
#       is.na(Xstart_full) ~ Inf,
#       .default = Xstart_full
#     )
#   ) %>%
#   rename(
#     `T`=tstop,
#     `V_early` = Xstart_early,
#     `V`=Xstart_full,
#     `J`= type
#   )
# 
# # Construct constraint matrix
# 
# Amat <- matrix(0, ncol=3, nrow=3)
# 
# for(i in 1:nrow(Amat)){
#   
#   if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#     
#     Amat[i,i]<--1
#     
#   }else{ # Force VE to wane (cannot increase over time)
#     
#     Amat[i,i] <- 1
#     
#   }
#   
# }
# 
# sieve_est_delta <- sieve_partially_linear_logistic(dat=delta_eligible, 
#                                                    V.early.name="V_early", 
#                                                    psi_delta = psi_d2_early, 
#                                                    monotone = TRUE, 
#                                                    Amat=Amat)
# 
# tmle_est_delta <- tmle_iterative(dat=delta_eligible, V.early.name="V_early", psi_delta=psi_d2_early, monotone=TRUE, Amat=Amat)
# 
# res_SLR <- data.frame(t=seq(0, max_t_delta/365.25, by=0.001) * 365.25)
# 
# grid <- cbind(as.numeric(seq(0,max_t_delta/365.25, by=0.001) < 14/365.25), 
#               as.numeric(seq(0, max_t_delta/365.25, by=0.001) >= 14/365.25),
#               ifelse(seq(0, max_t_delta/365.25, by=0.001) < 14/365.25, 0,
#                      seq(0, max_t_delta/365.25, by=0.001)))
# 
# res_SLR$VE_sieve <- 1-exp(grid %*% sieve_est_delta$beta)
# 
# CI_sieve <- monotone_CI_MC(beta = sieve_est_delta$beta_unconstr, 
#                vcov = sieve_est_delta$cov, 
#                Amat = sieve_est_delta$Amat, w = 1/(sieve_est_delta$se)^2,
#                grid = grid,
#                M=10000,
#                seed=47
#                )
# 
# res_SLR$VE_sieve_LCI <- 1-exp(CI_sieve[,2])
# res_SLR$VE_sieve_UCI <- 1-exp(CI_sieve[,1])
# 
# res_SLR$VE_TMLE <- 1-exp(grid %*% tmle_est_delta$beta)
# 
# CI_TMLE <- monotone_CI_MC(beta = tmle_est_delta$beta_unconstr, 
#                            vcov = tmle_est_delta$cov, 
#                            Amat = tmle_est_delta$Amat, w = 1/(tmle_est_delta$se)^2,
#                            grid = grid,
#                            M=10000,
#                            seed=47
# )
# 
# res_SLR$VE_TMLE_LCI <- 1-exp(CI_TMLE[,2])
# res_SLR$VE_TMLE_UCI <- 1-exp(CI_TMLE[,1])
# 
# res_SLR$type="Delta"
# 
# res_SLR <- res_SLR %>% 
#   pivot_longer(cols=c("VE_TMLE", "VE_TMLE_LCI", "VE_TMLE_UCI", "VE_sieve", "VE_sieve_LCI", "VE_sieve_UCI"), names_sep = "_", names_to=c("label", "method", "limit"), values_to="VE") %>%
#   mutate(limit = ifelse(is.na(limit), "est", limit)) %>%
#   dplyr::select(-label) %>%
#   rename(times=t)
# 
# res_SLR <- res_SLR %>% pivot_wider(names_from = "limit", values_from="VE")
# 
# tmp<-bind_rows(res_cox %>% mutate(method="Cox"), res_SLR %>% rename(VE = est, VE_lci = LCI, VE_uci=UCI))
# 
# ggplot()+
#   geom_ribbon(aes(x=times, ymin=VE_lci, ymax=VE_uci, fill=method), alpha=0.15, data=tmp %>% filter(type=="Delta"))+
#   geom_line(aes(x=times, y=VE, color=method), data=tmp %>% filter(type=="Delta"), linewidth=1.5)+
#   scale_color_manual(values = c("red", "purple", "blue"))+
#   scale_fill_manual(values = c("red", "purple", "blue"))+
#   scale_y_continuous(labels=scales::percent, limits=c(0,1))+
#   theme_minimal()
# 
# ########################
# # Fit Sieve & TMLE models to Omicron
# #########################
# 
# min_t_omicron <- as.numeric(difftime(dmy("14DEC2021"), dmy("23SEP2021")))
# max_t_omicron <- 194
# scale_f=365.25
# 
# omicron_eligible <- events %>% 
#   filter(tstop >= min_t_omicron & tstop <= max_t_omicron & type != "Delta") %>%
#   dplyr::select(Xstart_full, Xstart_early, tstop, type) %>%
#   mutate(tstop = tstop/scale_f, 
#          Xstart_early = Xstart_early/scale_f,
#          Xstart_full = Xstart_full/scale_f,
#          type = ifelse(type=="Omicron", 1, 0)) %>%
#   mutate(
#     Xstart_early = case_when(
#       is.na(Xstart_early) ~ Inf,
#       .default = Xstart_early
#     ),
#     Xstart_full = case_when(
#       is.na(Xstart_full) ~ Inf,
#       .default = Xstart_full
#     )
#   ) %>%
#   rename(
#     `T`=tstop,
#     `V_early` = Xstart_early,
#     `V`=Xstart_full,
#     `J`= type
#   )
# 
# max_tau = max(pmax(0, as.numeric(omicron_eligible$T-omicron_eligible$V)))*365.25
# 
# 
# #########
# # Try B-spline basis expansion
# ##########
# 
# Amat <- matrix(0, ncol=10, nrow=10)
# 
# for(i in 1:nrow(Amat)){
#   
#   if(i==1 | i==2){ # Force initial efficacy parameters to be negative
#     
#     Amat[i,i]<--1
#     
#   }else{ # Force VE to wane (cannot increase over time)
#     
#     Amat[i,(i-1)] <- -1
#     Amat[i,i] <- 1
#     
#   }
#   
# }
# 
# sieve_est_omicron <- sieve_partially_linear_logistic(dat=omicron_eligible, 
#                                                    V.early.name="V_early", 
#                                                    psi_delta = psi_bs_early, 
#                                                    monotone = TRUE, 
#                                                    Amat=Amat,
#                                                    df=8)
# 
# tmle_est_omicron <- tmle_iterative(dat=omicron_eligible, V.early.name="V_early", 
#                                    psi_delta=psi_bs_early, monotone=TRUE, Amat=Amat, smooth_r=TRUE, smooth_alpha=TRUE, df=8)
# 
# 
# grid <- cbind(as.numeric(seq(0, max_tau) <= 13),
#               #as.numeric(seq(0, max_tau) >= 14),
#                as.numeric(seq(0, max_tau) >= 14) * predict(sieve_est_omicron$basis, seq(0, max_tau)/365.25))
# 
# res_SLR <- data.frame(t=seq(0, max_tau))
# 
# res_SLR$VE_sieve <- 1-exp(grid %*% sieve_est_omicron$beta)
# 
# CI_sieve <- monotone_CI_MC(beta = sieve_est_omicron$beta_unconstr, 
#                            vcov = sieve_est_omicron$cov, 
#                            Amat = sieve_est_omicron$Amat, 
#                            w = 1/(sieve_est_omicron$se)^2,
#                            grid = grid,
#                            M=10000,
#                            seed=47
# )
# 
# res_SLR$VE_sieve_LCI <- 1-exp(CI_sieve[,2])
# res_SLR$VE_sieve_UCI <- 1-exp(CI_sieve[,1])
# 
# res_SLR$VE_TMLE <- 1-exp(grid %*% tmle_est_omicron$beta)
# 
# CI_TMLE <- monotone_CI_MC(beta = tmle_est_omicron$beta_unconstr, 
#                           vcov = tmle_est_omicron$cov, 
#                           Amat = tmle_est_omicron$Amat, w = 1/(tmle_est_omicron$se)^2,
#                           grid = grid,
#                           M=10000,
#                           seed=47
# )
# 
# res_SLR$VE_TMLE_LCI <- 1-exp(CI_TMLE[,2])
# res_SLR$VE_TMLE_UCI <- 1-exp(CI_TMLE[,1])
# 
# 
# # Construct constraint matrix
# 
# # Amat <- matrix(0, ncol=3, nrow=3)
# # 
# # for(i in 1:nrow(Amat)){
# #   
# #   if(i==1 | i==2){ # Force initial efficacy parameters to be negative
# #     
# #     Amat[i,i]<--1
# #     
# #   }else{ # Force VE to wane (cannot increase over time)
# #     
# #     Amat[i,i] <- 1
# #     
# #   }
# #   
# # }
# # 
# # sieve_est_omicron <- sieve_partially_linear_logistic(dat=omicron_eligible, 
# #                                                      V.early.name="V_early", 
# #                                                      psi_delta = psi_d2_early, monotone = TRUE, Amat=Amat)
# # 
# # tmle_est_omicron <- tmle_iterative(dat=omicron_eligible, V.early.name="V_early", 
# #                                    psi_delta=psi_d2_early, monotone=TRUE, Amat=Amat, smooth_r=TRUE, smooth_alpha=TRUE)
# # 
# # res_SLR <- data.frame(t=seq(0, max_tau/365.25, by=0.001) * 365.25)
# # 
# # grid <- cbind(as.numeric(seq(0,max_tau/365.25, by=0.001) < 14/365.25), 
# #               as.numeric(seq(0, max_tau/365.25, by=0.001) >= 14/365.25),
# #               ifelse(seq(0, max_tau/365.25, by=0.001) < 14/365.25, 0,
# #                      seq(0, max_tau/365.25, by=0.001)))
# # 
# # res_SLR$VE_sieve <- 1-exp(grid %*% sieve_est_omicron$beta)
# # 
# # CI_sieve <- monotone_CI_MC(beta = sieve_est_omicron$beta_unconstr, 
# #                            vcov = sieve_est_omicron$cov, 
# #                            Amat = sieve_est_omicron$Amat, w = 1/(sieve_est_omicron$se)^2,
# #                            grid = grid,
# #                            M=10000,
# #                            seed=47
# # )
# # 
# # res_SLR$VE_sieve_LCI <- 1-exp(CI_sieve[,2])
# # res_SLR$VE_sieve_UCI <- 1-exp(CI_sieve[,1])
# # 
# # res_SLR$VE_TMLE <- 1-exp(grid %*% tmle_est_omicron$beta)
# # 
# # CI_TMLE <- monotone_CI_MC(beta = tmle_est_omicron$beta_unconstr, 
# #                           vcov = tmle_est_omicron$cov, 
# #                           Amat = tmle_est_omicron$Amat, w = 1/(tmle_est_omicron$se)^2,
# #                           grid = grid,
# #                           M=10000,
# #                           seed=47
# # )
# # 
# # res_SLR$VE_TMLE_LCI <- 1-exp(CI_TMLE[,2])
# # res_SLR$VE_TMLE_UCI <- 1-exp(CI_TMLE[,1])
# 
# res_SLR$type="Omicron"
# 
# res_SLR <- res_SLR %>% 
#   pivot_longer(cols=c("VE_TMLE", "VE_TMLE_LCI", "VE_TMLE_UCI", "VE_sieve", "VE_sieve_LCI", "VE_sieve_UCI"), names_sep = "_", names_to=c("label", "method", "limit"), values_to="VE") %>%
#   mutate(limit = ifelse(is.na(limit), "est", limit)) %>%
#   dplyr::select(-label) %>%
#   rename(times=t)
# 
# res_SLR <- res_SLR %>% pivot_wider(names_from = "limit", values_from="VE")
# 
# tmp<-bind_rows(res_cox %>% mutate(method="Cox"), res_SLR %>% rename(VE = est, VE_lci = LCI, VE_uci=UCI))
# 
# ggplot()+
#   geom_ribbon(aes(x=times, ymin=pmax(0, VE_lci), ymax=VE_uci, fill=method), alpha=0.15, data=tmp %>% filter(type=="Omicron"))+
#   geom_line(aes(x=times, y=VE, color=method), data=tmp %>% filter(type=="Omicron"), linewidth=1.5)+
#   scale_color_manual(values = c("red", "purple", "blue"))+
#   scale_fill_manual(values = c("red", "purple", "blue"))+
#   scale_y_continuous(labels=scales::percent, limits=c(0,1), breaks=seq(0,1,by=0.1))+
#   labs(title="Omicron COVID-19", y = "VE (1-HR)", x="Days since boost")+
#   theme_minimal()+
#   theme(
#     #panel.grid=element_blank(),
#     axis.text=element_text(size=16),
#     axis.title=element_text(size=18),
#     legend.text=element_text(size=16),
#     legend.title=element_text(size=18),
#     legend.position="bottom",
#     strip.text=element_text(size=18)
#   )

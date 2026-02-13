#-------------------------------------------------------------------------------
# Name:        swash_test
# Purpose:     Tests and examples for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.3.0
# Last update: 2025-12-30 11:32
# Copyright (c) 2025 Thomas Wieland
#-------------------------------------------------------------------------------

library(lubridate)
library(sf)
library(spdep)

source("swash/R/swash.R")
# Loading swash code


# Switzerland:

load("swash/data/COVID19Cases_geoRegion.rda")

table(COVID19Cases_geoRegion$geoRegion)
table(COVID19Cases_geoRegion$datum)

COVID19Cases_geoRegion <-
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion %in% c("CH", "CHFL"),]
# CH = Switzerland all-over, CHFL = Switzerland and Liechtenstein all-over

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# First COVID-19 wave

COVID19Cases_geoRegion_balanced <- 
  is_balanced(
  data = COVID19Cases_geoRegion,
  col_cases = "entries",
  col_date = "datum",
  col_region = "geoRegion"
)
# Test whether "COVID19Cases_geoRegion" is balanced panel data 

COVID19Cases_geoRegion_balanced$data_balanced
# Balanced? TRUE or FALSE

CH_covidwave1 <- 
  swash (
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion"
  )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model

plot(CH_covidwave1)
# Plot of Swash-Backwash Model edges and total epidemic curve

plot_regions(
  CH_covidwave1,
  normalize_by_col = "pop",
  plot_rollmean = TRUE
  )


CH_covidwave1_confint <- 
  confint(
    CH_covidwave1, 
    iterations = 100
    )
# Bootstrap confidence intervals with 100 iterations

summary(CH_covidwave1_confint)
# Summary of confidence intervals

plot(CH_covidwave1_confint)
# Plot of confidence intervals

CH_covidwave1_growth <- growth(CH_covidwave1)
CH_covidwave1_growth
# Logistic growth models for sbm object CH_covidwave1

CH_covidwave1_initialgrowth_3weeks <- 
  growth_initial(
    CH_covidwave1,
    time_units = 21
    )
CH_covidwave1_initialgrowth_3weeks$results
# Exponential models for sbm object CH_covidwave1 
# initial growth in the first 3 weeks


# Austria:

load("swash/data/Oesterreich_Faelle.rda")

table(Oesterreich_Faelle$NUTS3)
table(Oesterreich_Faelle$Datum)

AT_covidwave1 <- 
  swash (
    data = Oesterreich_Faelle,
    col_cases = "Faelle",
    col_date = "Datum",
    col_region = "NUTS3"
  )
# Swash-Backwash Model for Austrian COVID19 cases
# Spatial aggregate: NUTS 3

summary(AT_covidwave1)

plot(AT_covidwave1)


AT_vs_CH <- 
  compare_countries(
    CH_covidwave1, 
    AT_covidwave1,
    country_names = c("Switzerland", "Austria"),
    iterations = 1000
    )

AT_vs_CH

plot(AT_vs_CH)


COVID19Cases_ZH <-
  COVID19Cases_geoRegion[
    (COVID19Cases_geoRegion$geoRegion == "ZH")
    & (COVID19Cases_geoRegion$sumTotal > 0),]
# COVID cases for Zurich

loggrowth_BS <- logistic_growth (
  y = COVID19Cases_ZH$sumTotal, 
  t = COVID19Cases_ZH$datum, 
  S = 3600,
  S_start = NULL, 
  S_end = NULL, 
  S_iterations = 10, 
  S_start_est_method = "bisect", 
  seq_by = 10,
  nls = TRUE
)
# Logistic growth model with stated saturation value

summary(loggrowth_BS)
# Summary of logistic growth model

plot(loggrowth_BS)
# Plot of logistic growth model

Rt_BS <- R_t(infections = COVID19Cases_ZH$entries)
Rt_BS
# Effective reproduction number


expgrowth_BS <- exponential_growth (
  y = COVID19Cases_ZH$sumTotal[1:28], 
  t = COVID19Cases_ZH$datum[1:28] 
)
# Exponential growth model for the first 4 weeks

summary(expgrowth_BS)
# Summary of exponential growth model

expgrowth_BS@doubling
# Doubling rate


load("swash/data/RKI_Corona_counties.rda")
# German counties (Source: Robert Koch Institute)

Corona_nbmat <- 
  nbmatrix (
    RKI_Corona_counties, 
    ID_col="AGS"
  )
# Creating neighborhood matrix

Corona_nbstat <- 
  nbstat (
    RKI_Corona_counties, 
    ID_col="AGS",
    link_data = RKI_Corona_counties, 
    data_ID_col = "AGS", 
    data_col = "EWZ", 
    func = "sum"
  )
Corona_nbstat$nbmat_data_aggregate
# Sum of population (EWZ) of neighboring counties


load("swash/data/did_fatalities_splm_coef.rda")
# Results of a difference-in-differences model

plot_coef_ci(
  point_estimates = did_fatalities_splm_coef$Estimate,
  confint_lower = did_fatalities_splm_coef$CI_lower_Bonferroni,
  confint_upper = did_fatalities_splm_coef$CI_upper_Bonferroni,
  coef_names = did_fatalities_splm_coef$Var,
  skipvars = c(
    "Alpha_share", 
    "lambda",
    "rho",
    "log(D_Infections_daily_7dsum_per100000_lag2weeks)",
    "vacc_cum_per100000_lag2weeks"
    ),
  lwd = 13,
  pch = 19,
  auto_color = TRUE
)
# Plot with point estimates and confidence intervals

plot_coef_ci(
  point_estimates = did_fatalities_splm_coef$Estimate,
  confint_lower = did_fatalities_splm_coef$CI_lower_Bonferroni,
  confint_upper = did_fatalities_splm_coef$CI_upper_Bonferroni,
  coef_names = did_fatalities_splm_coef$Var,
  p = did_fatalities_splm_coef$Pr_t_Bonferroni,
  skipvars = c(
    "Alpha_share", 
    "lambda",
    "rho",
    "log(D_Infections_daily_7dsum_per100000_lag2weeks)",
    "vacc_cum_per100000_lag2weeks"
  ),
  lwd = 13,
  pch = 19,
)
# Plot with point estimates and confidence intervals


load("swash/data/Infections.rda")
# Confirmed SARS-CoV-2 cases in Germany

plot_breakpoints(
  Infections, 
  log(infections_daily) ~ day,
  output.full = TRUE
  )
# Breakpoints for time series
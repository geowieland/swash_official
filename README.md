# swash: Health Geography Toolbox for Model-Based Analysis of Infections Panel Data

The R library `swash` provides a toolbox for quantitative analyses in health geography with a focus on the spatial spread of infectious diseases.
It bundles functions developed by the author between 2020 and 2023 for the quantitative analysis of (panel) infection data during the COVID-19 pandemic.
The aim was to consolidate these methods and analytical tools into a unified and coherent framework.
The target audience of this R package consists of researchers and practitioners in the fields of health geography, spatial epidemiology, and statistics.
Spread velocity may be analysed with the Swash-Backwash Model for the Single Epidemic Wave and corresponding functions for bootstrap confidence intervals, country comparison, and visualization of results.
Differences in epidemic growth between regions may be anaylsed using logistic growth models, exponential growth models, Hawkes processes and breakpoint analyses. All functionalities are accessed by the class `infpan` for infections panel data defined in this package, which is built from a `data.frame` provided by the user.


## Author

Thomas Wieland [ORCID](https://orcid.org/0000-0001-5168-9846) [EMail](mailto:geowieland@googlemail.com) 


## Availability

- 📦 CRAN: [swash](https://cran.r-project.org/package=swash)
- 💻 GitHub Repository: [swash_official](https://github.com/geowieland/swash_official)
- 📄 DOI (Zenodo): [10.5281/zenodo.18652150](https://doi.org/10.5281/zenodo.18652150)


## Citation

If you use this software, please cite:

Wieland, T. (2026). swash: Health Geography Toolbox for Model-Based Analysis of Infections Panel Data (Version 2.0.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.18652150


## Installation

To install the package from the CRAN repository, use:

```bash
install.packages("swash")
```

To install the package from GitHub: 

```bash
install.packages("remotes")
remotes::install_github("geowieland/swash_official")
```


## Features

The R library is a toolbox for quantitative analysis in health geography towards the spatial spread of infectious diseases.
In order to use all functionalities, the user should import her/his infections panel data using the function `load_infections_paneldata()`, 
which returns an instance of class `infpan`. The panel data is checked whether it is balanced and whether it includes missing values. 
From an `infpan` object, the user may utilize the following built-in analysis models and visualization functions:

- Swash-Backwash Model for the Single Epidemic Wave, including further analysis towards bootstrap-based inference and country comparison as well as visualization
- Growth Analysis with logistic and exponential growth models for cumulative or incremental infections, whereby the former is intended for the entire infection wave, and the latter for the initial phase of the infection wave; including visualization
- Hawkes process models for incremental infections; including visualization
- Breakpoints analysis using the Bai-Perron algorithm implemented in `strucchange::breakpoints`; including visualization
- Calculation of further epidemic indicators from the infections panel data such as the effective reproduction number
- Plots of infection curves by region

`infpan` objects and objects resulting from the functions mentioned above have `summary()` and `plot()` methods. 
All mentioned functions may be used stand-alone as well.


## Examples

```R
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_geoRegion <-
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion \%in\% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

infpan_CH <- load_infections_paneldata(
    data = COVID19Cases_geoRegion,
    col_cases = "entries",
    col_date = "datum",
    col_region = "geoRegion",
    other_cols = c(
      "Population" = "pop"
        ), 
    verbose = TRUE
  )
# Import as infections panel data set (class infpan)

is(infpan_CH)
# "infpan"

plot(
  infpan_CH,
  plot_rollmean = TRUE
  )
# Plot cases

infpan_CH <- calculate_Rt(
  infpan_CH,
  verbose = TRUE
  )
# Calculate effective reproduction number

infpan_CH <- calculate_rollmean(
  infpan_CH, 
  col_name = "RollingMean",
  verbose = TRUE
)
# Calculate rolling mean of cases as "RollingMean"

infpan_CH <- calculate_cum(
  infpan_CH, 
  col_name = "cumulatives",
  verbose = TRUE
)
# Calculate cumulative values of cases as "cumulatives"

infpan_CH <- calculate_incidence(
  infpan_CH, 
  col_name = "incidence",
  verbose = TRUE
)
# Calculate incidence of cases as "incidence"

summary(infpan_CH)
# Summary of infpan object

timestamps(infpan_CH)
# Time stamps of infpan object

CH_covidwave1_growth <- 
  growth(infpan_CH)
CH_covidwave1_growth
summary(CH_covidwave1_growth)
# Logistic growth models for infpan object infpan_CH

CH_covidwave1_initialgrowth_3weeks <- 
  growth_initial(
    infpan_CH,
    time_units = 21
  )
CH_covidwave1_initialgrowth_3weeks
summary(CH_covidwave1_initialgrowth_3weeks)
# Exponential models for infpan object CH_covidwave1 
# initial growth in the first 3 weeks

CH_covidwave1_Hawkes <- 
  growth_hawkes(infpan_CH)
CH_covidwave1_Hawkes
summary(CH_covidwave1_Hawkes)
# Hawkes process models for infpan object infpan_CH 

CH_covidwave1_breaks <- 
  growth_breaks(infpan_CH)
CH_covidwave1_breaks
summary(CH_covidwave1_breaks)
# Breakpoints for infpan object infpan_CH 

CH_covidwave1 <-
  swash(
    infpan_CH,
    verbose = TRUE
    )
# Swash-Backwash Model for Swiss COVID19 cases
# Spatial aggregate: NUTS 3 (cantons)

summary(CH_covidwave1)
# Summary of Swash-Backwash Model
```

See the /tests directory for usage examples of most of the included functions.


## Literature

Chowell, G., Simonsen, L., Viboud, C., & Yang, K. (2014). Is West Africa approaching a catastrophic phase or is the 2014 Ebola epidemic slowing down? Different models yield different answers for Liberia. *PLoS Currents*, 6. https://doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81

Chowell, G., Viboud, C., Hyman, J. M., & Simonsen, L. (2015). The Western Africa Ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. *PLOS Currents Outbreaks*. https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261

Cliff, A. D., & Haggett, P. (2006). A swash-backwash model of the single epidemic wave. *Journal of Geographical Systems*, 8(3), 227–252. https://doi.org/10.1007/s10109-006-0027-8

Li, M. Y. (2018). *An Introduction to Mathematical Modeling of Infectious Diseases*. Springer. https://doi.org/10.1007/978-3-319-72122-4

Nishiura, H., & Chowell, G. (2009). The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. In G. Chowell, J. M. Hyman, & L. M. A. Bettencourt (Eds.), *Mathematical and Statistical Estimation Approaches in Epidemiology* (pp. 103–121). Springer. https://doi.org/10.1007/978-90-481-2313-1_5

Pell, B., Kuang, Y., Viboud, C., & Chowell, G. (2018). Using phenomenological models for forecasting the 2015 Ebola challenge. *Epidemics*, 22, 62–70. https://doi.org/10.1016/j.epidem.2016.11.002

Rizoiu, M. A., Mishra, S., Kong, Q., Carman, M. & Xie, L. (2018). SIR-Hawkes: Linking Epidemic Models and Hawkes
Processes to Model Diffusions in Finite Populations. *Proceedings of the 2018 World Wide Web Conference. 
WWW’18*. Republic and Canton of Geneva, CHE: International World Wide Web Conferences
Steering Committee, p. 419–428. https://doi.org/10.1145/3178876.3186108

Smallman-Raynor, M. R., Cliff, A. D., & Stickler, P. J. (2022a). Meningococcal meningitis and coal mining in provincial England: Geographical perspectives on a major epidemic, 1929–33. *Geographical Analysis*, 54, 197–216. https://doi.org/10.1111/gean.12272

Smallman-Raynor, M. R., Cliff, A. D., & The COVID-19 Genomics UK (COG-UK) Consortium. (2022b). Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. *Epidemiology and Infection*, 150, e145. https://doi.org/10.1017/S0950268822001285

Viboud, C., Bjørnstad, O. N., Smith, D. L., Simonsen, L., Miller, M. A., & Grenfell, B. T. (2006). Synchrony, waves, and spatial hierarchies in the spread of influenza. *Science*, 312, 447–451. https://doi.org/10.1126/science.1125237

Wieland, T. (2020a). Flatten the curve! Modeling SARS-CoV-2/COVID-19 growth in Germany at the county level. *REGION*, 7(2), 43–83. https://doi.org/10.18335/region.v7i2.324

Wieland, T. (2020b). A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. *Safety Science*, 131, 104924. https://doi.org/10.1016/j.ssci.2020.104924

Wieland, T. (2022). Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. *European Journal of Geography*, 13(4), 18–33. https://doi.org/10.48088/ejg.t.wie.13.4.018.033

Wieland, T. (2025). Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: Results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. *Journal of Public Health: From Theory to Practice*, 33(11), 2497–2511. https://doi.org/10.1007/s10389-024-02218-x


## What's new (v2.0.2)

- General
  - Deprecation warnings with respect to version >=3.0.0

- Bugfixes
  - Updating URLs in help texts that are no longer valid
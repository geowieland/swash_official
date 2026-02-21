# swash: Implementation of the Swash-Backwash Model for the Single Epidemic Wave and additional functions in R

The `swash` package provides an implementation of the *Swash-Backwash Model for the Single Epidemic Wave* (Cliff and Haggett 2006) and additional functions for bootstrap confidence intervals, cross-country comparison, and data management. Other functions for spatio-temporal analysis and modeling of infectious diseases are included. The library is designed for researchers in health geography and similar fields.


## Author

Thomas Wieland [ORCID](https://orcid.org/0000-0001-5168-9846) [EMail](mailto:geowieland@googlemail.com) 


## Availability

- 📦 CRAN: [swash](https://cran.r-project.org/package=swash)
- 💻 GitHub Repository: [swash_official](https://github.com/geowieland/swash_official)
- 📄 DOI (Zenodo): [10.5281/zenodo.18652150](https://doi.org/10.5281/zenodo.18652150)


## Citation

If you use this software, please cite:

Wieland, T. (2026). swash: Implementation of the Swash-Backwash Model for the Single Epidemic Wave and additional functions in R (Version 1.3.3) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.18652150


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

- **Swash-Backwash Model**: 
  - Estimation of the Swash-Backwash Model (SBM) for panel data on regional infections
  - Bootstrap confidence intervals for SBM estimates
  - Cross-country comparison of model estimates (e.g., spatial reproduction numbers)
- **Spatio-temporal analysis of infections**:
  - Fitting logistic growth models for the whole epidemic wave
  - Fitting exponential growth models for the initial growth of the epidemic wave
  - Visual breakpoint analysis for time series
- **Other tools**:
  - Spatial statistics based of neighborhood matrix
  - Fit metrics for models with continuous and binary outcomes


## Examples

```R
data(COVID19Cases_geoRegion)
# Get SWISS COVID19 cases at NUTS 3 level

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[!COVID19Cases_geoRegion$geoRegion %in% c("CH", "CHFL"),]
# Exclude CH = Switzerland total and CHFL = Switzerland and Liechtenstein total

COVID19Cases_geoRegion <- 
  COVID19Cases_geoRegion[COVID19Cases_geoRegion$datum <= "2020-05-31",]
# Extract first COVID-19 wave

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
```

See the /tests directory for usage examples of most of the included functions.


## Literature

Cliff, A. D., & Haggett, P. (2006). A swash-backwash model of the single epidemic wave. *Journal of Geographical Systems*, 8(3), 227–252. https://doi.org/10.1007/s10109-006-0027-8

Chowell, G., Simonsen, L., Viboud, C., & Yang, K. (2014). Is West Africa approaching a catastrophic phase or is the 2014 Ebola epidemic slowing down? Different models yield different answers for Liberia. *PLoS Currents*, 6. https://doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81

Chowell, G., Viboud, C., Hyman, J. M., & Simonsen, L. (2015). The Western Africa Ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. *PLOS Currents Outbreaks*. https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261

Li, M. Y. (2018). *An Introduction to Mathematical Modeling of Infectious Diseases*. Springer. https://doi.org/10.1007/978-3-319-72122-4

Nishiura, H., & Chowell, G. (2009). The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. In G. Chowell, J. M. Hyman, & L. M. A. Bettencourt (Eds.), *Mathematical and Statistical Estimation Approaches in Epidemiology* (pp. 103–121). Springer. https://doi.org/10.1007/978-90-481-2313-1_5

Pell, B., Kuang, Y., Viboud, C., & Chowell, G. (2018). Using phenomenological models for forecasting the 2015 Ebola challenge. *Epidemics*, 22, 62–70. https://doi.org/10.1016/j.epidem.2016.11.002

Smallman-Raynor, M. R., Cliff, A. D., & Stickler, P. J. (2022a). Meningococcal meningitis and coal mining in provincial England: Geographical perspectives on a major epidemic, 1929–33. *Geographical Analysis*, 54, 197–216. https://doi.org/10.1111/gean.12272

Smallman-Raynor, M. R., Cliff, A. D., & The COVID-19 Genomics UK (COG-UK) Consortium. (2022b). Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. *Epidemiology and Infection*, 150, e145. https://doi.org/10.1017/S0950268822001285

Viboud, C., Bjørnstad, O. N., Smith, D. L., Simonsen, L., Miller, M. A., & Grenfell, B. T. (2006). Synchrony, waves, and spatial hierarchies in the spread of influenza. *Science*, 312, 447–451. https://doi.org/10.1126/science.1125237

Wieland, T. (2020). Flatten the curve! Modeling SARS-CoV-2/COVID-19 growth in Germany at the county level. *REGION*, 7(2), 43–83. https://doi.org/10.18335/region.v7i2.324

Wieland, T. (2020). A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. *Safety Science*, 131, 104924. https://doi.org/10.1016/j.ssci.2020.104924

Wieland, T. (2022). Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. *European Journal of Geography*, 13(4), 18–33. https://doi.org/10.48088/ejg.t.wie.13.4.018.033

Wieland, T. (2025). Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: Results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. *Journal of Public Health: From Theory to Practice*, 33(11), 2497–2511. https://doi.org/10.1007/s10389-024-02218-x


## What's new (v1.3.3)
- New features:
  - Fit metrics for logistic and exponential growth models.
- Bugfixes:
  - Bug in calculation in metrics() fixed
  - Check for vectors lengths in logistic_growth() and exponential_growth()
  - Deprecation warnings with respect to version >=2.0.0.
# Infections data

## Definitions

Let $$N$$ be the number of observational units (e.g., individuals, firms, regions), and $$T$$ the number of time points (e.g., days, weeks, years) a variable is observed. In the present context, the observed variable is typically [new or cumulative infections](infections-indicators.md#new-infections-and-cumulative-infections), and the observational units are countries or regions within a country. Based on this, we make the following simplified distinction between three types of data[1]:

### Cross-sectional data

If different observational units ($$N>1$$) are examined at a single point in time ($$T=1$$), then it is a case of *cross-sectional data*. Such a data set has $$N$$ rows.

### Time-series data

*Time-series data* means that a single observational unit ($$N=1$$) is observed over multiple time periods ($$T>1$$), which means that the data set has $$T$$ rows.

### Panel data

*Panel data* are a combination of cross-sectional and time-series da: they consist of data in which the same units of analysis ($$N > 1$$) are observed repeatedly ($$T > 1$$), typically at regular time intervals. (Complete) panel data sets have $$N \cdot T$$ rows.

When the same observational units are examined at two points in time ($$N>1$$, $$T=2$$), with the first point occurring before an event or intervention and the second point after it, this is referred to as *pre-post data*. When there are more than two points in time for data collection, this may be referred to as *multi-period panel data*[2][3].

If a panel dataset is complete - which means that there is an observation for every combination of observational unit and time point - it is referred to as a *balanced panel* with $$N \cdot T$$ rows. If even a single observation is missing for any observational unit or time point, it is an *unbalanced panel* including fewer than $$N \cdot T$$ rows[1]. Missing values with respect to specific time points in time-series or panel data may be filled in via imputation[4].

## Infections data in spatial analyses of epidemic spreads

The [Swash-Backwash Model of the Single Epidemic Wave](swash-backwash.md) requires multi-period panel data on new infections. It does not strictly require the panel data to be balanced in order for it to be calculated. However, the results will be biased if individual data points are missing. 

Phenomenological epidemic models, such as the [logistic growth model](logistic-growth.md), for one observational unit and a given number of time points may be calculated using incomplete data as well; however, obtaining an unbiased estimate requires a complete time series. The same applies if, for example, growth rates are to be compared across regions, necessitating multi-period panel data. 

Spatial statistics, e.g., to identify significant clusters, are calculated only for cross-sectional data. They are therefore not affected by unbalanced panels, however, a problem may occur if specific data points are missing for the particular moment in time under consideration.

## Examples of public infections data

During the COVID-19 pandemic, for example, the Swiss Federal Office of Public Health (FOPH) published panel data on daily and cumulative infections, broken down by canton ([NUTS 3](https://en.wikipedia.org/wiki/Nomenclature_of_Territorial_Units_for_Statistics))[3][5]. Many countries published COVID-19 data in this or a similar way. In Germany, data with a higher level of resolution were published, namely a dataset containing daily case reports with further information[6], which may be aggregated into panel data at the level of counties (NUTS 3) and days[3][7]. There were also scientific projects that collected COVID-19 case data from various countries and aggregated and consolidated it at the national or regional level, for example, the data from Johns Hopkins University[8] or the *COVID-19 European regional tracker*[9].

The German Robert Koch Institute collects data on cases of infectious diseases and the detection of pathogens that are subject to mandatory reporting under the German Infection Protection Act. These cases are aggregated anonymously and can be accessed via the *SurvStat* service. The data can be aggregated according to various criteria, including administrative regions (e.g., counties) and time periods (in this instance, calendar weeks), allowing users to generate panel data for specific infectious diseases or pathogens[10][11].

## References

[1] Stock JH, Watson MW (2019) 
*Introduction to Econometrics*. Fourth Edition. 
New York: Pearson.

[2] Nianogo RA, Benmarhnia T, O'Neill S (2023)
A comparison of quasi-experimental methods with data before and after an intervention: an introduction for epidemiologists and a simulation study. 
*International Journal of Epidemiology* 52(5), 1522–1533. [10.1093/ije/dyad032](https://doi.org/10.1093/ije/dyad032) 

[3] Wieland T (2025) 
Assessing the effectiveness of non-pharmaceutical interventions in the SARS-CoV-2 pandemic: Results of a natural experiment regarding Baden-Württemberg (Germany) and Switzerland in the second infection wave. 
*Journal of Public Health: From Theory to Practice*, 33(11), 2497–2511. [10.1007/s10389-024-02218-x](https://doi.org/10.1007/s10389-024-02218-x)

[4] Jahangiri M, Kazemnejad A, Goldfeld KS, Daneshpour MS, Mostafaei S, Khalili D, Moghadas MR, Akbarzadeh M. (2023)
A wide range of missing imputation approaches in longitudinal data: a simulation study and real data analysis. 
*BMC Medical Research Methodology* 23(1), 161. [10.1186/s12874-023-01968-8](https://doi.org/10.1186/s12874-023-01968-8)

[5] Bundesamt für Gesundheit (2023) 
COVID-19 dashboard source data – DailyIncomingData Schema. 
https://covid19.admin.ch#/definitions/DailyIncomingData (Accessed June 28, 2023)

[6] Robert Koch Institut (2020) 
*Tabelle mit den aktuellen Covid-19 Infektionen pro Tag (Zeitreihe)*.  
https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/-dd4580c810204019a7b8eb3e0b329dd60/data (Accessed June 28, 2020)

[7] Wieland T (2020) 
Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
*REGION* 7(2), 43–83. [10.18335/region.v7i2.324](https://doi.org/10.18335/region.v7i2.324)

[8] Johns Hopkins University (2020) 
New Cases of COVID-19 In World Countries.
https://coronavirus.jhu.edu/data/new-cases (accessed May 11, 2020)

[9] Naqvi A (2021) 
COVID-19 European regional tracker. 
*Scientific Data* 8, 181. [10.1038/s41597-021-00950-7](https://doi.org/10.1038/s41597-021-00950-7)

[10] Faensen D, Krause G. (2004) 
SurvStat@RKI – a web-based solution to query surveillance data in Germany. 
*Eurosurveillance* 8(22),2477. [10.2807/esw.08.22.02477-en](https://doi.org/10.2807/esw.08.22.02477-en)

[11] Robert Koch Institut (2026) 
SURVSTAT@RKI 2.0. 
https://survstat.rki.de/Content/Query/Create.aspx (accessed October 24, 2026)
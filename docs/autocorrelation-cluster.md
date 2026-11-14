# Statistics for spatial autocorrelation and clusters

## Background

Spatial dependence means that the value of a variable at a specific location is determined (in part) by the values ​​of the same or other variables at neighboring locations, due to the existence of a (statistical and causal) relationship. Spatial dependence is measured by statistical indicators of spatial autocorrelation[1]. In the context of epidemic spreads, the mathematical determination of spatial dependence is primarily used to investigate clusters and hotspots of a specific infectious disease (e.g., in terms of [incidence](infections-indicators.md#incidence), deaths per capita, or [excess mortality](https://en.wikipedia.org/wiki/Excess_mortality)). Communicable diseases often exhibit spatial dependence: since spatial interactions between subregions of a larger area decrease as distance increases ([distance decay](https://en.wikipedia.org/wiki/Distance_decay)), clusters with high or relatively low case densities frequently emerge. Commuter flows play a significant role in the transmission of many infectious diseases. Furthermore, clusters also arise from individual super-spreading events[2][3][4][5][6]. 

## Spatial autocorrelation

Spatial autocorrelation refers to the degree to which observations of a variable at nearby geographical locations are statistically related to each other, i.e., it is the correlation with the values ​​of neighboring or nearby regions. It describes the spatial pattern of a variable, indicating whether similar values tend to cluster together (positive spatial autocorrelation), dissimilar values are located near each other (negative spatial autocorrelation), or spatial patterns occur purely randomly (no spatial autocorrelation). Spatial autocorrelation is an important concept in [spatial analysis](https://en.wikipedia.org/wiki/Spatial_analysis), which is typically used to identifying underlying geographical patterns and dependencies in the data[1]. For the analysis of spatial autocorrelation or the detection of clusters there are both global and local approaches[1][2].

## Global statistics

Global spatial autocorrelation statistics yield a single indicator value representing the global autocorrelation across all subareas. The best-known ones are *Moran's I* and *Getis-Ord G*, as shown in the following[7][8]:

### Moran's I

Moran's I measures global spatial autocorrelation by determining whether similar values are spatially clustered, randomly distributed, or dispersed across a study area:

$$I = \frac{n}{W} \cdot \frac{\sum_i \sum_j w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_i (x_i - \bar{x})^2}$$

where $$I$$ is the Moran's I statistic, $$n$$ is the number of spatial units, $$w_{ij}$$ is the [spatial weight between units $$i$$ and $$j$$](neighborhood-matrix.md#spatial-weights-matrix), $$W$$ is the sum of all spatial weights, $$x_i$$ and $$x_j$$ are the observed values for spatial units $$i$$ and $$j$$, respectively, and $$\bar{x}$$ is the mean of all observed values.

The expected value under the null hypothesis (= there is no spatial autocorrelation) equals:

$$E(I) = \frac{-1}{n-1}$$

where $$E(I)$$ is the expected value of $$I$$ under the null hypothesis, and $$n$$ is the number of spatial units within the analysis.

The variance of $$I$$ is calculated as follows:

$$Var(I) = E(I^2) - [E(I)]^2$$

where $$E(I)$$ is the expected value of Moran's I under the null hypothesis of spatial randomness, and $$E(I^2)$$ is the expected squared value of Moran's I.

For the significance test, the value of $$I$$ is converted into the following test statistic:

$$Z = \frac{I - E(I)}{\sqrt{Var(I)}}$$

where $$I$$ is the observed Moran's I, $$E(I)$$ is its expected value under the null hypothesis, and $$Var(I)$$ is the variance of Moran's I. The resulting Z-score is used to assess whether the observed spatial autocorrelation is statistically significant.

### Getis-Ord G

Getis-Ord G identifies global spatial concentrations of high or low values, indicating whether hotspots or cold spots exist across the study area:

$$G = \frac{\sum_i \sum_j w_{ij} x_i x_j}{\sum_i \sum_j x_i x_j}$$

where $$G$$ is the Getis-Ord statistic, $$w_{ij}$$ is the [spatial weight between units $$i$$ and $$j$$](neighborhood-matrix.md#spatial-weights-matrix), and $$x_i$$ and $$x_j$$ are the observed attribute values for spatial units $$i$$ and $$j$$, respectively.

The expected value under the null hypothesis (= there is no spatial autocorrelation) equals:

$$E(G) = \frac{W}{n(n-1)}$$

where $$W$$ is the sum of all spatial weights and $$n$$ is the number of spatial units.

The variance of $$G$$ is calculated as follows:

$$Var(G) = E(G^2) - [E(G)]^2$$

where $$E(G^2)$$ is the expected squared value of the Getis-Ord statistic under the null hypothesis of spatial randomness.

For the significance test, the value of $$G$$ is converted into the following test statistic:

$$Z = \frac{G - E(G)}{\sqrt{Var(G)}}$$

where $$G$$ is the observed Getis-Ord statistic, $$E(G)$$ is its expectd value under the null hypothesis, and $$Var(G)$$ is its variance. The resulting Z-score indicates whether high or low values are significantly spatially clustered. A positive Z-score indicates clustering of high values, whereas a negative Z-score indicates clustering of low values.

## Local statistics

The latter approaches calculate indicators for any spatial object. These are typically called *LISA (Local Indicators of Spatial Association)*[9][10]:

### Local Getis-Ord Gi Statistic $$G_i^*$$

Getis-Ord $$G_i^*$$ identifies local spatial clusters of high or low values by detecting statistically significant hotspots and cold spots at specific locations within a study area:

$$G_i^* = \frac{\sum_j w_{ij} x_j}{\sum_j x_j}$$

where $$G_i^*$$ is the local $$G$$ statistic for spatial unit $$i$$, $$w_{ij}$$ is the [spatial weight between spatial units $$i$$ and $$j$$](neighborhood-matrix.md#spatial-weights-matrix), and $$x_j$$ is the observed attribute value for spatial unit $$j$$. The statistic measures the local spatial concentration of high or low values around a specific spatial unit $$i$$, *including the unit $$i$$ itself*.

The expected value under the null hypothesis (= there is no clustering surrounding region $$i$$) equals:

$$E(G_i^*) = \frac{\sum_j w_{ij}}{n}$$

where $$n$$ is the total number of spatial units and $$w_{ij}$$ represents the spatial weights.

The variance of $$G_i^*$$ is calculated as follows:

$$Var(G_i^*) = E[(G_i^*)^2] - [E(G_i^*)]^2$$

where $$E[(G_i^*)^2]$$ is the expected squared value of the local Getis-Ord Gi* statistic under the null hypothesis of spatial randomness.

For the significance test, the value of $$G_i^*$$ is converted into the following test statistic:

$$Z(G_i^*) = \frac{G_i^* - E(G_i^*)}{\sqrt{Var(G_i^*)}}$$

where $$G_i^*$$ is the observed local Getis-Ord statistic, $$E(G_i^*)$$ is its expected value under the null hypothesis, and $$Var(G_i^*)$$ is its variance. The resulting $$Z$$ score indicates whether high or low values are significantly spatially clustered around location $$i$$. A significantly positive $$Z$$ score indicates a *hotspot*, whereas a significantly negative $$Z$$ score indicates a *cold spot*.

## Further notes

In the context of the SARS-CoV-2/COVID-19 outbreak that began in 2020, several studies examined the spatial distribution of SARS-CoV-2 infections and/or COVID-19-associated deaths, employing indicators of spatial autocorrelation such as [Moran's I](#morans-i). In all cases, these analyses revealed a greater or lesser degree of spatial dependence, with several clearly identifiable clusters[6][11][12][13]. In one instance, an analysis of German counties during the first SARS-CoV-2 wave revealed that both COVID-19-associated deaths per capita and the intrinsic growth rates derived from [logistic growth models](logistic-growth.md) exhibit spatial autocorrelation[14].


## References

[1] Griffith D (2009) 
Spatial Autocorrelation. 
In: Kitchin R, Thrift N (eds), *International Encyclopedia of Human Geography*, 308-316.

[2] Lessler J, Salje H, Grabowski MK, Cummings DAT (2016) 
Measuring Spatial Dependence for Infectious Disease Epidemiology.
*PLoS ONE* 11(5), e0155249. [10.1371/journal.pone.0155249](https://doi.org/10.1371/journal.pone.0155249)

[3] Rotejanaprasert C, Chinpong K, Lawson AB, Chienwichai P, Maude RJ (2024)
Evaluation and comparison of spatial cluster detection methods for improved decision making of disease surveillance: 
a case study of national dengue surveillance in Thailand. 
*BMC Medical Research Methodology* 24, 14. [10.1186/s12874-023-02135-9](https://doi.org/10.1186/s12874-023-02135-9)

[4] Charaudeau S, Pakdaman K, Böelle PY (2014). 
Commuter Mobility and the Spread of Infectious Diseases: Application to Influenza in France. 
*PLOS ONE* 9(1), 1-9. [10.1371/journal.pone.0083002](https://doi.org/10.1371/journal.pone.0083002)

[5] Charu V, Zeger S, Gog J, Bjørnstad ON, Kissler S, Simonsen L, Grenfell BT, Viboud, C (2017). 
Human mobility and the spatial transmission of influenza in the United States. 
*PLOS Computational Biology*, 13, 1-23. https://doi.org/10.1371/journal.pcbi.1005382.

[6] Wieland T (2022) 
Spatial patterns of excess mortality in the first year of the COVID-19 pandemic in Germany. 
*European Journal of Geography*, 13(4), 18–33. [10.48088/ejg.t.wie.13.4.018.033](https://doi.org/10.48088/ejg.t.wie.13.4.018.033)

[7] Getis A, Ord JK (1992)
The analysis of spatial association by use of distance statistics. 
*Geographical Analysis* 24 (3), 189–206. [10.1111/j.1538-4632.1992.tb00261.x](https://doi.org/10.1111/j.1538-4632.1992.tb00261.x)

[8] Moran PAP (1950) 
Notes on Continuous Stochastic Phenomena. 
*Biometrika* 37(1/2), 17-23. [10.2307/2332142](https://doi.org/10.2307/2332142)

[9] Ord JK, Getis A (1995) 
Local Spatial Autocorrelation Statistics: Distributional Issues and an Application. 
*Geographical Analysis* 27, 286-306. [10.1111/j.1538-4632.1995.tb00912.x](https://doi.org/10.1111/j.1538-4632.1995.tb00912.x)

[10] Danek T, Weglinska E, Zareba M (2022) 
The influence of meteorological factors and terrain on air pollution concentration and migration: a geostatistical case study from Krakow, Poland. 
*Scientific Reports* 12, 11050. [10.1038/s41598-022-15160-3](https://doi.org/10.1038/s41598-022-15160-3)

[11] Bourdin S, Jeanne L, Nadou F, Noiret G (2021)
Does lockdown work? A spatial analysis of the spread and concentration of Covid-19 in Italy.
*Regional Studies* 55(7), 1182–1193. [10.1080/00343404.2021.1887471](https://doi.org/10.1080/00343404.2021.1887471)

[12] Saffary T, Adegboye OA, Gayawan E, Elfaki F, Kuddus MA, Saffary R (2020) 
Analysis of COVID-19 Cases' Spatial Dependence in US Counties Reveals Health Inequalities.
*Frontiers in Public Health* 8, 579190. [10.3389/fpubh.2020.579190](https://doi.org/10.3389/fpubh.2020.579190)

[13] Wang Q, Dong W, Yang K, Ren Z, Huang D, Zhang P, Wang J (2021)
Temporal and spatial analysis of COVID-19 transmission in China and its influencing factors.
*International Journal of Infectious Diseases* 105, 675-685. [10.1016/j.ijid.2021.03.014](https://doi.org/10.1016/j.ijid.2021.03.014)

[14] Wieland T (2020) 
Flatten the curve! Modeling SARS-CoV-2/COVID-19 growth in Germany at the county level. 
*REGION*, 7(2), 43–83. [10.18335/region.v7i2.324](https://doi.org/10.18335/region.v7i2.324)
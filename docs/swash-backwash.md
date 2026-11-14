# Swash-Backwash Model of the Single Epidemic Wave

## Background

The *Swash-Backwash Model of the Single Epidemic Wave* (in short: SBM) was developed by Cliff and Haggett (2006)[1]. It is the spatial equivalent of the classic epidemiological [Susceptible-Infected-Recovered (SIR) Model](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)#The_SIR_model), which divides the population into the not yet infected (susceptible), the infected, and the recovered, and models the development of these three groups over the course of a wave of an infectious disease. The model draws its terminology from coastal geomorphology, where waves advance up a beach (*swash*) and subsequently retreat (*backwash*). 

In simple terms, the SBM analyses the spread of an infectious disease across a country's regions, distinguishing between regions where no cases have yet occurred (*susceptible areas*), those with active infections (*infected areas*), and those where cases are no longer occurring (*recovered areas*). The transmission speed between regions is represented by the *spatial reproduction number* $$R_{0A}$$, which refers to the [basic reproduction number $$R_0$$](exponential-growth.md#basic-reproduction-number-and-doubling-time) in epidemiology[1][2][3]. The analysis of a spatial epidemic spread using the SBM requires [infections panel data](infections-data.md#panel-data) of [new infections](infections-indicators.md#new-infections-and-cumulative-infections) (e.g., daily, weekly or monthly).
  

## Concept and mathematical formulation

We have an area $$A$$ (e.g., country) which consists of $$J$$ ($$j = 1,2,...,J$$) subregions (e.g., counties, federal states). For these regions, infection panel data are available for $$T$$ ($$t = 1,2,...,T$$) discrete time periods. We denote $$N_{jt}$$ as the number of new infection cases in region $$j$$ at time $$t$$, and the total number of infection cases during the epidemic wave across all regions and time points is $$N = \sum_{j=1}^J \sum_{t=1}^T N_{jt}$$. 

Over the course of an infection wave, the regions within the total area of investigation fall into one of three mutually exclusive states: regions where no infections have (yet) occurred at time $$t$$ (*susceptible areas*), regions where infections are occurring at time $$t$$ (*infected areas*), and regions where infections are no longer occurring at time $$t$$ (*recovered areas*). Just as the sum of the susceptible, infected, and recovered populations in the SIR model equals the total population, the sum of the three area integrals corresponds to the entire area[1]:

$$S_A + I_A + R_A = 1 = A$$

where $$S_A$$ contains all susceptible regions, $$I_A$$ contains all infected regions, and $$R_A$$ contains all recovered regions within the total area of investigation, $$A$$.

For a given region $$j$$ in which an infection first occurred (or was first detected), we refer to the point in time $$t$$ when this took place (e.g., day, week, month) as the *leading edge* (LE). The point in time when the last infection was reported is called the *following edge* (FE). For both edges, a time-weighted arithmetic mean is defined, which indicates the average time of arrival or, respectively, the average time of decay of the epidemic wave across all regions[1][2]:

$$\bar{t}_{LE} = \frac{1}{J_{LE}} \sum_{t=1}^T t \cdot N_{t,LE}$$

$$\bar{t}_{FE} = \frac{1}{J_{FE}} \sum_{t=1}^T t \cdot N_{t,FE}$$

where $$\bar{t}_{LE}$$ and $$\bar{t}_{LE}$$ are time-weighted arithmetic means, which represent the average time of arrival and departure of the epidemic wave across the regions, and $$N_{t,LE}$$ and $$N_{t,FE}$$ is the number of regions whose leading edge and following edge, respectively, occurred in time period $$t$$.

An initial measure of the rate of spatial spread is the difference between the two:

$$\Delta t = \bar{t}_{FE}-\bar{t}_{LE}$$

The susceptible (proportion of the study area at risk of infection), infected (proportion of the study area which is infected), and recovered integrals can be calculated from the rate of spread as follows[1][2]:

$$S_A = \frac{\bar{t}_{LE}-1}{T}$$

$$I_A = \frac{\bar{t}_{FE}}{T} - S_A$$

$$R_A = 1 - (S_A + I_A)$$

where $$S_A$$ is the susceptible integral, $$I_A$$ is the infected integral, and $$R_A$$ is the recovered integral.

Drawing on considerations regarding the mathematical formulation of the basic reproduction number, Cliff and Haggett derive an analogous spatial reproduction number, which they define as follows: "In spatial terms, $$A$$, the spatial reproduction number,
$$R_{0A}$$, is the average number of secondary districts produced from one infected district in a virgin area"[1]. It is calculated as follows[1][2]:

$$R_{0A} = \frac{1-S_A}{1-R_A}$$

where $$R_{0A}$$ is the spatial reproduction number.

## Bootstrap inference

*Bootstrap sampling* is a *resampling* technique in which multiple new samples of the same size are repeatedly drawn with replacement from an original dataset. It allows the statistical properties of a given estimator to be empirically assessed. Bootstrap sampling forms the basis of *bootstrap inference*, which includes generating [confidence intervals](https://en.wikipedia.org/wiki/Confidence_interval) without making strong assumptions about the underlying population distribution ([nonparametric method](https://en.wikipedia.org/wiki/Nonparametric_statistics)). It is especially valuable when the true distribution is unknown or when traditional statistical methods are not applicable. Confidence intervals may be obtained from the quantiles of the bootstrap distribution. For a 95% confidence interval, the lower and upper bounds correspond to the 2.5th and 97.5th percentiles, respectively[4][5][6]. This is the case for the results of the Swash-Backwash Model (e.g., spatial reproduction number $$R_{0A}$$), as its basic form does not include any statistical inference (see above).

Bootstrap inference for the SBM results may be generated as follows:

1. Define a [significance level](https://en.wikipedia.org/wiki/Statistical_significance) $$\alpha$$, e.g., $$\alpha = 0.05$$ 
2. Define a number of $$B$$ bootstrap samples and the sample size $$n$$, more precisely: the proportion $$s$$ of of the $$J$$ regions which is to be drawn from the original dataset (e.g., 0.8), which leads to $$n = s \cdot J \cdot T$$ (entity-based bootstrap sampling)
3. Draw $$B$$ samples  (with replacement), each consisting of $$s \cdot J$$ regions, from the original data and, for each of these, calculate the SBM (see above)
4. For the SBM parameter(s) (e.g., $$R_{0A}$$), calculate the $$\frac{\alpha}{2}$$ and $$1-\frac{\alpha}{2}$$ percentiles from the bootstrap distribution, which represent the $$1 - \alpha \cdot 100$$ confidence intervals

## Nonparametric test for differences between two countries

The SBM does not offer a way to perform inferential statistical comparisons of SBM parameter results between different areas (e.g., countries) or between infectious diseases within the same country. For comparing the means of two samples of approximately normally distributed variables, several mean-comparison methods based on the t-distribution (["t-tests"](https://en.wikipedia.org/wiki/Student%27s_t-test)) are available. However, these are not applicable to estimators with unknown distributions. In such cases, nonparametric statistical methods can also be used to examine whether there is a significant difference in a parameter $$\hat{\theta}$$ between two independent groups. While this parameter may be a mean value, the method is, in principle, applicable to any parameter whatsoever. The key idea behind this is to calculate the differences in the parameter for each bootstrap sample and then examine the distribution of these difference values[4][7][8]. 

Suppose we have two (meaningfully comparable) countries, $$A$$ and $$B$$, experiencing an epidemic spread. These countries include $$J_A$$ and $$J_B$$ regions, respectively. To examine whether there is a significant difference in a specific parameter $$\hat{\theta}$$ between these countries (e.g., $$R_{0A}$$ from the SBM), one can proceed as follows[7][8]:

1. Define a [significance level](https://en.wikipedia.org/wiki/Statistical_significance) $$\alpha$$, e.g., $$\alpha = 0.05$$ 
2. Define a number of $$B$$ bootstrap samples ($$b = 1,2,...,B$$) and the proportion $$s$$ of of the $$J_A$$ and $$J_B$$ regions which is to be drawn from the original dataset (e.g., 0.8), which leads to $$n_A = s \cdot J_A \cdot T$$ and $$n_B = s \cdot J_B \cdot T$$
3. Draw $$B$$ samples (with replacement) from countries $$A$$ and $$B$$, respectively, each consisting of $$s \cdot J_A$$ and $$s \cdot J_B$$ regions, respectively, from the original data and, for each of these, calculate the SBM (see above)
4. Combine the two arrays $$\hat{\theta}^{*}_{A,b}$$ and $$\hat{\theta}^{*}_{B,b}$$ and calculate the difference between the estimators: $$D_b^* = \hat{\theta}^{*}_{A,b} - \hat{\theta}^{*}_{B,b}$$
5. Calculate the $$\frac{\alpha}{2}$$ and $$1 - \frac{\alpha}{2}$$ percentiles from the distribution of $$D_b^*$$, which represent the $$1 - \alpha \cdot 100$$ confidence intervals of $$D_b^*$$  

If both percentiles are below or above zero, we may interpret the difference in parameter $$\hat{\theta}^{*}$$ between the countries as significant at the level of $$\alpha$$.

## Further notes

Several historical outbreaks have already been analyzed using this model[1][2]. The model was also used for a comparative analysis of the spatial spread rates of different SARS-CoV-2 variants[3].

## References

[1] Cliff, A. D. & Haggett, P. (2006). 
A swash-backwash model of the single epidemic wave. 
*Journal of Geographical Systems*, 8(3), 227–252. 
[10.1007/s10109-006-0027-8](https://doi.org/10.1007/s10109-006-0027-8)

[2] Smallman-Raynor, M. R., Cliff, A. D., & Stickler, P. J. (2022). 
Meningococcal meningitis and coal mining in provincial England: Geographical perspectives on a major epidemic, 1929–33. 
*Geographical Analysis*, 54, 197–216. [10.1111/gean.12272](https://doi.org/10.1111/gean.12272)

[3] Smallman-Raynor, M. R., Cliff, A. D., & The COVID-19 Genomics UK (COG-UK) Consortium. (2022). 
Spatial growth rate of emerging SARS-CoV-2 lineages in England, September 2020–December 2021. 
*Epidemiology and Infection*, 150, e145. [10.1017/S0950268822001285](https://doi.org/10.1017/S0950268822001285)

[4] Dixon PM (2006) Bootstrap Resampling. In *Encyclopedia of Environmetrics* (eds A.H. El-Shaarawi and W.W. Piegorsch). [10.1002/9780470057339.vab028](https://doi.org/10.1002/9780470057339.vab028)

[5] Efron B, Tibshirani R (1986)
Bootstrap Methods for Standard Errors, Confidence Intervals, and Other Measures of Statistical Accuracy.
*Statistical Science* 1(1), 54-75. [10.1214/ss/1177013815](https://doi.org/10.1214/ss/1177013815) 

[6] Rousselet GA, Pernet CR, Wilcox RR (2021) 
The Percentile Bootstrap: A Primer With Step-by-Step Instructions in R. 
*Advances in Methods and Practices in Psychological Science* 4(1), 1-10. [10.1177/2515245920911881](https://doi.org/10.1177/2515245920911881)

[7] Hesterberg TC (2015) 
What Teachers Should Know About the Bootstrap: Resampling in the Undergraduate Statistics Curriculum. 
*The American Statistician* 69(4), 371–386. [10.1080/00031305.2015.1089789](https://doi.org/10.1080/00031305.2015.1089789)

[8] Çetinkaya-Rundel M, Hardin J (2026)
*Introduction to Modern Statistics*. Second Edition.
https://openintrostat.github.io/ims/
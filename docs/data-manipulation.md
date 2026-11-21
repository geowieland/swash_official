# Data manipulation

## Background

When dealing with time series or panel data, there may occur some difficulties such as missing values, non-positive values, or strong fluctuations. In the case of infection data, this can be due to, for example, reporting delays or very high (e.g., daily) measurement frequency. 

## Imputation of missing values in time series

Missing values ​​can be filled in using linear imputation. In this process, the missing values are estimated using linear interpolation between the nearest known observations before and after the missing value. For a missing value between two known values $$y_{t_0}$$ and $$y_{t_1}$$, the imputed value is calculated as[<a href="#ref-1">1</a>][<a href="#ref-2">2</a>]:

$$\hat{y}_t = y_{t_0} + \frac{t-t_0}{t_1-t_0} (y_{t_1}-y_{t_0})$$

where $$y_{t_0}$$ and $$y_{t_1}$$ are the time points of the known observations, and $$\hat{y}_t$$ is the imputed value. This method assumes a linear change between the two observed values.

## Handling zeros with log1p transformation 

Non-positive values cannot always be processed in the same way as positive values. The occurrence of negative values in infections time series would indicate erroneous data (as negative infections are impossible). However, values equal to zero are realistic, for example, in cases without any (reported) infection. These values ​​pose a problem when the time series needs to be log-transformed (e.g., for the OLS estimation of an intrinsically linear model such as the [exponential growth model](exponential-growth.md)).

The *log1p transformation* applies a logarithm after adding a constant to each data point[<a href="#ref-3">3</a>][<a href="#ref-4">4</a>]:

$$y' = \log(y+c)$$

Although $$c=1$$ is commonly used, other constants may be chosen depending on the data. The condition $$y+c>0$$ must hold for all values. Please note that this method introduces a bias into the data (and into any estimates derived from them) the magnitude of which depends on the number of zero values ​​and the range of the variable being transformed[<a href="#ref-3">3</a>][<a href="#ref-4">4</a>].

## Smoothing of time series using rolling averages or sums

Rolling (or: moving) sums and means are used to smooth short-term fluctuations in time series, such as daily time series, making underlying trends and temporal patterns easier to identify. Rolling sums and rolling means are calculated over a moving window of consecutive observations in a time series. For a time point $$t$$, they are defined as[<a href="#ref-2">2</a>]:

$$S_t = \sum_{i=0}^{k-1} y_{t-i}$$

and

$$\bar{y}_t = \frac{1}{k}\sum_{i=0}^{k-1} y_{t-i}$$

where $$y_{t-i}$$ denotes the observation at time $$t-1$$, $$k$$ is the window size, $$S_t$$ is the rolling sum, and $$\bar{y}_t$$ is the rolling mean. The rolling window can be left-aligned, centered, or right-aligned, depending on how the window is positioned relative to time point. The formulas above correspond to a right-aligned window, which includes the current observation and the preceding observations.

For example, during the SARS-CoV-2/COVID-19 pandemic, the German federal public health institute (Robert Koch Institut, RKI) calculated the [incidence](infections-indicators.md#incidence) of reported SARS-CoV-2 infection cases as a rolling sum of cases over the last 7 days per 100,000 inhabitants[<a href="#ref-5">5</a>]:

$$Inc7_{t} = \frac{\sum_t^{t-6} N_{t}}{pop_{t}} \cdot 100000$$

where $$Inc7_{t}$$ is the 7-day-incidence at time point $$t$$, $$N_{t}$$ is the number of new reported infection cases at time $$t$$, and $$pop_{t}$$ is the at time point $$t$$.


## References

<a id="ref-1"></a>
[1] Jahangiri M, Kazemnejad A, Goldfeld KS, Daneshpour MS, Mostafaei S, Khalili D, Moghadas MR, Akbarzadeh M. (2023)
A wide range of missing imputation approaches in longitudinal data: a simulation study and real data analysis. 
*BMC Medical Research Methodology* 23(1), 161. [10.1186/s12874-023-01968-8](https://doi.org/10.1186/s12874-023-01968-8)

<a id="ref-2"></a>
[2] Zeileis A, Grothendieck G (2005) 
zoo: S3 Infrastructure for Regular and Irregular Time Series.
*Journal of Statistical Software* 14(6), 1–27.
[10.18637/jss.v014.i06](https://doi.org/10.18637/jss.v014.i06)

<a id="ref-3"></a>
[3] Mullahy J, Norton EC (2024)
Why Transform Y? The Pitfalls of Transformed Regressions with a Mass at Zero.
*Oxford Bulletin of Economics and Statistics* 86(2), 417-447. [10.1111/obes.12583](https://doi.org/10.1111/obes.12583)

<a id="ref-4"></a>
[4] O'Hara RB, Kotze DJ (2010)
Do not log-transform count data. 
*Methods in Ecology and Evolution* 1, 118-122. [10.1111/j.2041-210X.2010.00021.x](https://doi.org/10.1111/j.2041-210X.2010.00021.x)

<a id="ref-5"></a>
[5] Frank C, Hellenbrand W, Sievers C (2021)
Betrachtung der reiseassoziierten COVID-19-Fälle im Sommer 2020 unter Berücksichtigung der Schulferien, Reisetätigkeit und Testkapazitäten.
*Epidemiologisches Bulletin* 8, 10-23. [10.25646/7955.3](https://doi.org/10.25646/7955.3)
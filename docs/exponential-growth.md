# Exponential growth model

## Background

Just like the [logistic growth model](logistic-growth.md), the exponential growth model is one of the *phenomenological models* for epidemics. In contrast to *mechanistic models* that aim to directly represent mechanisms of infection and spread (e.g., the [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)#The_SIR_model) and its extensions), phenomenological models of infection waves are statistical models that depict the growth of new infections (e.g., daily or weekly) or cumulative infections over time. A phenomenological model for epidemic modeling has the advantage of being less "data-hungry" and requiring no assumptions about the epidemiological parameters of the infectious disease whose wave is being studied[<a href="#ref-1">1</a>][<a href="#ref-2">2</a>][<a href="#ref-3">3</a>]. 

An exponential growth model is characterized by a constant relative growth rate and a constant *doubling time*. This means that the number of new infections increases by the same percentage over equal time intervals. As a result, not only does the total number of infections rise, but the absolute increase in new infections also grows from one period to the next[<a href="#ref-2">2</a>][<a href="#ref-4">4</a>]. The empirical parameterization of a single exponential growth model requires a [time series](infections-data.md#time-series-data) of [new infections](infections-indicators.md#new-infections-and-cumulative-infections) for the area under study; consequently, comparing exponential growth rates across multiple countries or regions requires [panel data](infections-data.md#panel-data) on new infections.

## Mathematical formulation

The exponential growth function is formalized as[<a href="#ref-4">4</a>][<a href="#ref-5">5</a>]:

$$N_t = N_0 \cdot e^{r \cdot t}$$

where $$N_t$$ is the number of [new infections](infections-indicators.md#new-infections-and-cumulative-infections) at time point $$t$$, $$t$$ is the time counter, and $$N_0$$ and $$r$$ are empirical parameters, with $$r$$ being the exponential growth rate per time unit (e.g., day, week) and $$N_0$$ being the start value (infections at time point 0).

## Basic reproduction number and doubling time

The *basic reproduction number* is intended to represent the number of secondary infections resulting from a primarily infected individual at the onset of an epidemic spread, assuming a completely susceptible population. The *doubling time* indicates the period of time within which the number of infection cases doubles[<a href="#ref-4">4</a>][<a href="#ref-5">5</a>]. Both may be calculated based on the exponential growth rate[<a href="#ref-5">5</a>][<a href="#ref-6">6</a>]:

$$R_0 = e^{r \cdot T_g}$$

and

$$T_d = \frac{\ln{2}}{r}$$

where $$R_0$$ is the basic reproduction number, $$T_d$$ is the doubling time, $$T_g$$ is the [*generation time*](https://en.wikipedia.org/wiki/Generation_time), and $$r$$ is the exponential growth rate ($$r>0$$) per time unit (e.g., day, week).

## Empirical parametrization

The exponential model may be estimated using [Nonlinear Least Squares](https://en.wikipedia.org/wiki/Non-linear_least_squares), given empirical data on new infections (e.g., per day or per week). As an *intrinsically linear model*, the exponential model may be linearized and estimated using [Ordinary Least Squares](https://en.wikipedia.org/wiki/Ordinary_least_squares)[<a href="#ref-6">6</a>][<a href="#ref-7">7</a>]:

$$\ln{N_t} = \ln{\hat{N}_0} + \hat{r} \cdot t + \ln{\epsilon_t}$$

where $$\hat{N}_0$$ and $$\hat{r}$$ are the empirical estimates of $$N_0$$ and $$r$$, and $$\ln{\epsilon_t}$$ is the additive stochastic disturbance term with respect to time point $$t$$.

This linearization requires that every individual empirical value $$N_t$$ is greater than 0. If values ​​of 0 occur in a time series of infections, all values may be [increased by a constant](data-manipulation.md#handling-zeros-with-log1p-transformation) (e.g., 1), though this will alter the estimate. Values ​​less than 0 (i.e., negative infections) are not possible, provided the input data are correct.

## Further notes

The exponential growth model for infections is used to analyze the *initial* phase of an epidemic, specifically to calculate rates of spread at the outset. In epidemiological literature, this is sometimes referred to as an "exponential phase" of an infection wave. The exponential growth rate and the basic reproduction number $$R_0$$ derived from it are regarded as rough estimates and are intended for the situation of a completely susceptible population[<a href="#ref-2">2</a>][<a href="#ref-4">4</a>]. Some studies also find sub-exponential growth even during the initial phase of an outbreak, and therefore classify the exponential growth model as less suitable[<a href="#ref-8">8</a>].

The exponential growth model is also utilized to identify differences in regional growth rates within a country or an area which is exposed to an epidemic spread[<a href="#ref-8">8</a>][<a href="#ref-9">9</a>][<a href="#ref-10">10</a>].

In the given context, the generation time describes the time span between the moment a person becomes infected and the moment that person infects another individual. It is frequently specified as a distribution[<a href="#ref-11">11</a>].

$$R_0$$ and $$T_g$$ can only be meaningfully calculated when $$r > 0$$. However, if $$r < 0$$ (exponential decay), then the doubling time corresponds to the *halving time* (or *half-life*), when it is calculated using $$|r|$$ instead of $$r$$[<a href="#ref-12">12</a>].

## References

<a id="ref-1"></a>
[1] Chowell G, Simonsen L, Viboud C, Yang K (2014) 
Is West Africa Approaching a Catastrophic Phase or is the 2014 Ebola Epidemic Slowing Down? Different Models Yield Different Answers for Liberia. 
*PLoS currents* 6. [10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81](https://dx.doi.org/10.1371/currents.outbreaks.b4690859d91684da963dc40e00f3da81)

<a id="ref-2"></a>
[2] Chowell G, Sattenspiel L, Bansal S, Viboud C (2016) 
Mathematical models to characterize early epidemic growth: A review. 
*Physics of Life Reviews* 18, 66-97. [10.1016/j.plrev.2016.07.005](https://doi.org/10.1016/j.plrev.2016.07.005)

<a id="ref-3"></a>
[3] Wieland T (2020) 
Flatten the Curve! Modeling SARS-CoV-2/COVID-19 Growth in Germany at the County Level. 
*REGION* 7(2), 43–83. [10.18335/region.v7i2.324](https://doi.org/10.18335/region.v7i2.324)

<a id="ref-4"></a>
[4] Nishiura H, Chowell G (2009) 
The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. 
In G Chowell, JM Hyman, LMA Bettencourt (Eds.), *Mathematical and Statistical Estimation Approaches in Epidemiology*, 
103–121. Springer. [10.1007/978-90-481-2313-1_5](https://doi.org/10.1007/978-90-481-2313-1_5)

<a id="ref-5"></a>
[5] Anzai A, Nishiura H (2022)
Doubling time of infectious diseases.
*Journal of Theoretical Biology* 554, 111278. [10.1016/j.jtbi.2022.111278](https://doi.org/10.1016/j.jtbi.2022.111278)

<a id="ref-6"></a>
[6] Bonifazi G, Lista L, Menasce D, Mezzetto M, Pedrini D, Spighi R, Zoccoli A (2021) 
A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. 
*The European Physical Journal Plus* 136, 386. 
[10.1140/epjp/s13360-021-01339-6](https://doi.org/10.1140/epjp/s13360-021-01339-6)

<a id="ref-7"></a>
[7] Rawlings JO, Pantula SG, Dickey DA (1998)
*Applied Regression Analysis: A Research Tool. Second Edition*.
New York: Springer.
[10.1007/b98890](https://doi.org/10.1007/b98890)

<a id="ref-8"></a>
[8] Chowell G, Viboud C, Hyman JM, Simonsen L (2015) 
The Western Africa Ebola virus disease epidemic exhibits both global exponential and local polynomial growth rates. 
*PLOS Currents Outbreaks*, ecurrents.outbreaks.8b55f4bad99ac5c5db3663e916803261. [10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261](https://doi.org/10.1371/currents.outbreaks.8b55f4bad99ac5c5db3663e916803261)

<a id="ref-9"></a>
[9] Messner W, Payson SE (2020)
Variation in COVID-19 outbreaks at the US state and county levels.
*Public Health* 187, 15-18.
[10.1016/j.puhe.2020.07.035](https://doi.org/10.1016/j.puhe.2020.07.035)

<a id="ref-10"></a>
[10] Burghardt K, Guo S, Lerman K (2022)
Unequal impact and spatial aggregation distort COVID-19 growth rates. 
*Philosophical Transactions A* 380 (2214), 20210122. [10.1098/rsta.2021.0122](https://doi.org/10.1098/rsta.2021.0122)

<a id="ref-11"></a>
[11] Svensson A (2007)
A note on generation times in epidemic models.
*Mathematical Biosciences* 208(1), 300-311. [10.1016/j.mbs.2006.10.010](https://doi.org/10.1016/j.mbs.2006.10.010)

<a id="ref-12"></a>
[12] Hunt R (1990) *Basic Growth Analysis*. Dordrecht: Springer. https://doi.org/10.1007/978-94-010-9117-6_3
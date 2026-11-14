# Detecting breakpoints in time series

## Background

In the simplest case, the growth of infections is represented as a function of time, for example, in an [exponential growth model](exponential-growth.md). This requires [time series](infections-data.md#time-series-data) of [new or cumulative infections](infections-indicators.md#new-infections-and-cumulative-infections). To determine turning points - that is, time points at which the growth trajectory changes - we can, as an alternative to the [logistic growth model](logistic-growth.md), identify structural changes (or breaks) within the time series. The Bai-Perron algorithm[1] may be used to identify one or more structural breaks within a time series model.

## Structural changes in time series

Consider a linear regression model for a time series with $$T$$ time points or observations:

$$y_t = \alpha + \beta \cdot t + \epsilon_t$$

where $$y_t$$ is the dependent variable, $$t$$ is the time index $$t = 1,2,...,T$$, $$\alpha$$ and $$\beta$$ are the regression coefficients to be estimated (intercept and slope, respectively), and $$epsilon_t$$ denotes the stochastic error term.

A structural change in the time series occurs when, from a given time point $$t_C$$ onwards, one or more regression coefficients (intercept and/or slope) differ significantly from their values before $$t_C$$. To test the hypothesis that the regression coefficients differ significantly before and after a specific, a priori defined time point $$t_C$$, one can use the *Chow test*. The Chow test examines the null hypothesis that the regression coefficients are identical across the two subperiods against the alternative hypothesis that at least one coefficient differs. To represent a structural break at time $$t_C$$, the time-series model is extended to include a dummy variable that indicates time points from $$t_C$$ onwards, as well as an interaction variable between this dummy and $$t$$[2][3].

## The Bai-Perron algorithm

The Bai-Perron algorithm functions differently from the Chow test and also serves a different purpose: In the absence of a hypothesis, it may be used to identify one or more structural breaks within a time-series model as described above[1][4]. 

The detection and dating of structural breaks within a time series includes following the subsequent steps[1][4][5]:

1. The regression model is specified based on the full time series, assuming constant regression coefficients over the entire observation period. This model serves as the reference specification against which potential structural changes are evaluated.
2. Potential structural breaks are identified using statistical tests, including fluctuation tests (recursive estimation tests) and F-statistics. These tests assess whether the regression coefficients remain stable over time or whether allowing for structural changes significantly improves the model fit.
3. For a given number of structural breaks $$M$$, the Bai-Perron algorithm determines the optimal breakpoint locations by minimizing the residual sum of squares (RSS) across all possible segmentations, subject to a minimum segment length. This procedure divides the time series into $$M+1$$ segments, within which the regression coefficients are assumed to remain constant.
4. To determine the optimal number of breakpoints, alternative model specifications containing zero to $$M_{test}$$ structural breaks (corresponding to one to $$M_{test}+1$$ segments) are estimated and compared. The optimal number of breaks is selected based on model selection criteria, for example, the [Bayesian Information Criterion (BIC)](https://en.wikipedia.org/wiki/Bayesian_information_criterion), with lower values indicating a preferred model.
5. Finally, the estimated breakpoint locations are assigned to their corresponding observations in the time series. Confidence intervals may be calculated to account for the statistical uncertainty surrounding the estimated breakpoint dates.

Consequently, the original time series is divided into $$S$$ distinct segments ($$s = 1,2,...,S$$), with $$S = M+1$$, shown here with a linear regression model. The regression coefficients are assumed to remain constant within each segment but are allowed to differ across segments[5]:

$$y_t = \alpha_s + \beta_s \cdot t + \epsilon_t$$

where $$\alpha_s$$ and $$\beta_s$$ are the regression coefficents of model segment $$s$$. 

If, for example, in segment $$s =1$$, the slope $$\beta_1$$ is positive and in segment $$s = 2$$, the slope $$\beta_2$$ is negative, then $$y_t$$ increases in segment 1 and declines in segment 2.

For the sake of simplicity, the examples presented above use a linear regression model with untransformed variables (which is unrelated to the algorithm used for detecting structural changes). However, when infections are the dependent variable, it is advisable to use an intrinsically linear model (e.g., an [exponential growth model](exponential-growth.md))[5][8]:

$$\ln{N_t} = \ln{\hat{N}_{0_s}} + \hat{r}_s \cdot t + \ln{\epsilon_t}$$

where $$N_t$$ is the number of [new infections](infections-indicators.md#new-infections-and-cumulative-infections) at time point $$t$$, and $$\hat{N}_{0_s}$$ and $$\hat{r}_s$$ are the regression coefficents of model segment $$s$$.

## Further notes

In the context of the SARS-CoV-2/COVID-19 spread, several studies investigated the effectiveness of lockdowns by identifying such structural breaks in infection time series in order to check whether they coincide with the imposement of lockdowns. These studies used either the Bai-Perron algorithm[5][6][7] or similar algorithms for the detection of breakpoints[8][9][10][11].

## References

[1] Bai J, Perron P (2003) 
Computation and analysis of multiple structural change models. 
*Journal of Applied Econometrics* 18(1), 1-22. [10.1002/jae.659](https://doi.org/10.1002/jae.659)

[2] Chow GC (1960) 
Tests of Equality Between Sets of Coefficients in Two Linear Regressions. 
*Econometrica* 28(3), 591–605. [10.2307/1910133](https://doi.org/10.2307/1910133)

[3] Lee FY (1974) 
On the Dummy Variable Technique and Covariance Analysis in Testing Equality Among Sets of Coefficients in Linear Regressions: An Expository Note. 
*The Journal of Financial and Quantitative Analysis* 9(3), 491–495. [10.2307/2329876](https://doi.org/10.2307/2329876)

[4] Zeileis C, Kleiber W, Krämer K, Hornik, K (2003) 
Testing and dating of structural changes in practice. 
*Computational Statistics & Data Analysis* 44(1-2), 109-123. [10.1016/S0167-9473(03)00030-6](https://doi.org/10.1016/S0167-9473(03)00030-6)

[5] Wieland T (2020)
A phenomenological approach to assessing the effectiveness of COVID-19 related nonpharmaceutical interventions in Germany. 
*Safety Science* 131, 104924. [10.1016/j.ssci.2020.104924](https://doi.org/10.1016/j.ssci.2020.104924)

[6] Hartl T, Wälde K, Weber E (2020)
Measuring the impact of the German public shutdown on the spread of COVID-19. 
*VoxEU*, April 14 2020. https://cepr.org/voxeu/columns/measuring-impact-german-public-shutdown-spread-covid-19

[7] Wälde K, Weber E, Hartl T (2020)
Measuring the impact of the German public shutdown on the spread of COVID-19. 
*Covid Economics* 1, 25–32. https://cepr.org/publications/covid-economics-issue-1#392514_392877_390183

[8] Küchenhoff H, Günther F, Höhle M, Bender A (2021) 
Analysis of the early COVID-19 epidemic curve in Germany by regression models with change points.
*Epidemiology and Infection* 149, e68, 1–7. [10.1017/S0950268821000558](https://doi.org/10.1017/S0950268821000558)

[9] Boudou M, O'hAiseadha C, Garvey P, O'Dwyer J, Hynds P (2021) 
Breakpoint modelling of temporal associations between non-pharmaceutical interventions and symptomatic COVID-19 incidence in the Republic of Ireland.
*PLoS ONE* 16(7), e0255254. [10.1371/journal.pone.0255254](https://doi.org/10.1371/journal.pone.0255254)

[10] Santamaría L, Hortal J (2020)
Chasing the ghost of infection past: identifying thresholds of change during the COVID-19 infection in Spain. 
*Epidemiology and Infection* 148, e282, 1–12. [10.1017/S0950268820002782](https://doi.org/10.1017/S0950268820002782)

[11] Santamaría L, Hortal J (2021)
COVID-19 effective reproduction number dropped during Spain's nationwide dropdown, then spiked at lower-incidence regions.
*Science of The Total Environment* 751, 142257. [10.1016/j.scitotenv.2020.142257](https://doi.org/10.1016/j.scitotenv.2020.142257)
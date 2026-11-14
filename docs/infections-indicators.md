# Indicators of infections over time

## New infections and cumulative infections

New infections are the number of individuals who become infected during a specified time interval. Cumulative infections are the total number of individuals who have been infected up to a given point in time, including all new infections recorded in previous intervals[1]. In this handbook, we consistently use the following notations:

- New infections at time $$t$$ (e.g., days, weeks) are denoted as $$N_t$$
- Cumulative infections at time $$t$$ (e.g., days, weeks) are denoted as $$C_t$$

## Incidence

The term "incidence" is not always used consistently. A unifying definition is as follows: incidence is the number of infections or cases of disease within a defined population and over a defined time period, which can be expressed either as an absolute number or as a density value (e.g., per capita)[1]. Here, we use incidence as a density value, e.g.:  

$$Inc_{t} = \frac{N_{t}}{pop_{t}} \cdot k$$

where $$Inc_{t}$$ is the incidence of a disease at time point $$t$$, $$N_{t}$$ is the number of new infections (or, more general, cases) at time point $$t$$, $$pop_{t}$$ is the population at time point $$t$$, and $$k$$ is a scaling factor.

For example, during the SARS-CoV-2/COVID-19 pandemic, the German federal public health institute (Robert Koch Institut, RKI) calculated the incidence of reported SARS-CoV-2 infection cases as a [rolling sum](data-manipulation.md#smoothing-of-time-series-using-rolling-averages-or-sums) of cases over the last 7 days per 100,000 inhabitants[2]:

$$Inc7_{t} = \frac{\sum_t^{t-6} N_{t}}{pop_{t}} \cdot 100000$$

where $$Inc7_{t}$$ is the 7-day-incidence at time point $$t$$, $$N_{t}$$ is the number of new reported infection cases at time $$t$$, and $$pop_{t}$$ is the at time point $$t$$.

## Effective reproduction number $$R_t$$

The effective reproduction number describes the average number of secondary infections originating from an infected individual at time point $$t$$. A value greater than one indicates an increase in infections, while a value equal to one corresponds to linear growth, and a value less than one corresponds to a decline in infections[3]. Here, we use the effective reproduction number — analogous to the approach taken by Germany's federal public health institute (Robert Koch Institut, RKI) — as an empirical metric calculated directly from the new infections[4][5]:

$$R_t = \frac{\sum_{j=0}^{k-1} N_{t-j}} {\sum_{j=k}^{2k-1} N_{t-j}}$$

where $$R_t$$ is the effective reproduction number at time point (day) $$t$$, $$N_t$$ Number of new infections on day $$t$$, $$k$$ is the length of the time window in days, and $$j$$ is the summation index used to iterate over the days within each time window.

## Further notes

The effective reproduction number has a conceptual link to the [basic reproduction number](exponential-growth.md#basic-reproduction-number-and-doubling-time), although the latter is not intended as an indicator which changes over time.


## References

[1] Porta M (2008) 
*A Dictionary of Epidemiology*. Fifth Edition. 
Oxford: Oxford University Press.

[2] Frank C, Hellenbrand W, Sievers C (2021)
Betrachtung der reiseassoziierten COVID-19-Fälle im Sommer 2020 unter Berücksichtigung der Schulferien, Reisetätigkeit und Testkapazitäten.
*Epidemiologisches Bulletin* 8, 10-23. [10.25646/7955.3](https://doi.org/10.25646/7955.3)

[3] Nishiura H, Chowell G (2009) 
The effective reproduction number as a prelude to statistical estimation of time-dependent epidemic trends. 
In G Chowell, JM Hyman, LMA Bettencourt (Eds.), *Mathematical and Statistical Estimation Approaches in Epidemiology*, 
103–121. Springer. [10.1007/978-90-481-2313-1_5](https://doi.org/10.1007/978-90-481-2313-1_5)

[4] Bonifazi G, Lista L, Menasce D, Mezzetto M, Pedrini D, Spighi R, Zoccoli A (2021) 
A simplified estimate of the effective reproduction number Rt using its relation with the doubling time and application to Italian COVID-19 data. 
*The European Physical Journal Plus* 136, 386. 
[10.1140/epjp/s13360-021-01339-6](https://doi.org/10.1140/epjp/s13360-021-01339-6)

[5] an der Heiden M, Hamouda O (2020) 
Schätzung der aktuellen Entwicklung der SARS-CoV-2-Epidemie in Deutschland - Nowcasting.
*Epidemiologisches Bulletin* 17, 10-15. [10.25646/6692](https://doi.org/10.25646/6692)
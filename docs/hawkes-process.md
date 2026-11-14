# Hawkes process

## Background

The Hawkes process is named after its developer[1] and has many applications, including in financial market econometrics, geophysics, and epidemiology[2]. It is a *self-exciting* point process that can be understood as an extension of the classical Poisson process, in which the occurrence of an event temporarily increases the intensity, and thus the likelihood, of subsequent events. This self-exciting effect gradually diminishes over time (*decay*)[1][2][3][4]. 

In the context of an epidemic spread, Hawkes processes are useful for modeling the intensity of infectious diseases by distinguishing between the contributions of exogenous and endogenous processes, that is, between a *background transmission rate* and local secondary infections. The relative weight of these two components is disease-specific and may vary significantly depending on the particular infectious disease[3][4]. While the [exponential](exponential-growth.md) and [logistic growth models](logistic-growth.md) are classified as *phenomenological models* for infection growth (in contrast to mechanistic models such as the [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)#The_SIR_model) and its extensions), the Hawkes process is referred to as a *semi-mechanistic model*[3]. The empirical parameterization of a single Hawkes process model requires a [time series](infections-data.md#time-series-data) of [new infections](infections-indicators.md#new-infections-and-cumulative-infections) for the area under study; consequently, comparing Hawkes process parameters across multiple countries or regions requires [panel data](infections-data.md#panel-data) on new infections.

## Mathematical formulation

The conditional intensity function of a Hawkes process with an exponential kernel is represented mathematically as follows[1][4]:

$$\lambda(t_i) = \mu + \alpha \cdot \sum_{j=1}^{i-1} y_j \cdot e^{\left(-\beta(t_i-t_j)\right)}$$

where $$\lambda(t_i)$$ denotes the conditional intensity of the event under study at time point $$t_i$$, $$y_j$$ is the number of events under study occurring at time point $$j$$ ($$j < t$$), $$t_i$$ and $$t_j$$ are time points $$t$$ and $$j$$, respectively, $$\mu$$ is the background intensity, $$\alpha$$ is the excitation parameter, and $$\beta$$ is the decay parameter. 

When the Hawkes process is applied to the growth of new infections, the parameters are to be interpreted as follows[3][4][5]: 
- *Background rate* $$\mu$$: equals the contribution to the intensity from *exogenous* sources, which means that it is independent from the self-exciting process, e.g., external importations of the infectious disease.
- *Excitation* $$\alpha$$: Determines the strength of the self-exciting effect. Each infection occurring within the system contributes $$\alpha$$ to the initial increase in the infection intensity, all other things being equal. This parameter therefore captures the endogenous component of infection growth, which reflects person-to-person transmission within the region under study. This contribution subsequently decreases exponentially over time.
- *Decay* $$\beta$$: Determines the rate at which the influence of previous infections on the infection intensity decreases exponentially over time. A higher value of $$\beta$$ indicates a faster decay of the self-exciting influence, whereas a lower value indicates a more persistent influence of previous infections.

From the conditional intensity equation with an exponential kernel, we may calculate the *branching ratio*[2][6]:

$$\gamma = \frac{\alpha}{\beta}$$

where $$gamma$$ is the branching ratio.

## Empirical parametrization

Unlike, for example, the exponential growth model, the Hawkes equation shown above is not an intrinsically linear model, which means that it cannot be linearized via a transformation and parameterized using Ordinary Least Squares (OLS). For this reason, a [Maximum Likelihood Estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) is used in most applications[2][3][5][6].

## Further notes

It has been shown in the epidemiological literature that the Hawkes process and the SIR model can be transformed into one another[3][5].

The illustration above shows the basic Hawkes process with an exponential kernel. There are many other variants of this[2][4].


## References

[1] Hawkes AG (1971)
Spectra of some self-exciting and mutually exciting point processes.
*Biometrika* 58(1), 83-90. [10.2307/2334319](https://doi.org/10.2307/2334319)

[2] Shlomovich L, Cohen EAK, Adams N, Patel L (2022)
Parameter Estimation of Binned Hawkes Processes.
*Journal of Computational and Graphical Statistics* 31(4), 990-1000. [10.1080/10618600.2022.2050247](https://doi.org/10.1080/10618600.2022.2050247)

[3] Unwin HJT, Routledge I, Flaxman S, Rizoiu M-A, Lai S, Cohen J, Weiss DJ, Mishra S, Bhatt S (2021) 
Using Hawkes Processes to model imported and local malaria cases in near-elimination settings. 
*PLoS Computational Biology* 17(4), e1008830. [10.1371/journal.pcbi.1008830](https://doi.org/10.1371/journal.pcbi.1008830)

[4] Laub PJ, Taimre T, Pollett PK (2015)
Hawkes Processes.
https://people.smp.uq.edu.au/PhilipPollett/papers/Hawkes.pdf

[5] Rizoiu MA, Mishra S, Kong Q, Carman M, Xie L (2018)
SIR-Hawkes: Linking Epidemic Models and Hawkes Processes to Model Diffusions in Finite Populations. 
*Proceedings of the 2018 World Wide Web Conference. WWW’18*, p. 419–428. Republic and Canton of Geneva, CHE: International World Wide Web Conferences Steering Committee. [10.1145/3178876.3186108](https://doi.org/10.1145/3178876.3186108)

[6] Hardiman SJ, Bouchaud J (2014) 
Branching-ratio approximation for the self-exciting Hawkes process.
*Physical Review E* 90, 062807. [10.1103/PhysRevE.90.062807](https://doi.org/10.1103/PhysRevE.90.062807)
---
bibliography: references.bib
---

# How reliable is the estimated phenology ?

## Introduction

We developped a method, based on the study of the productivity (proportion of juveniles among all caught individuals [@cuchot]), which allow the estimation of : - the breeding timing (xmid parameter, proxy for phenology) - the breeding success (asymptote parameter, final proportion of juvenile in the sampled population) - the variance in breeding phenology (scale parameter)

The aim of this project is to investigate how well estimated is the phenology (using sigmoid bayesian approach) when parameters associated with the phenology (mean and variance of laying dates, survival, breeding success, number of breeding attempts ...) vary. For this, we want to simulate breeding phenology of a (bird) population, sample it, according to a specified sampling design, and then compare the estimated values with the simulated ones.

(Define space, time)

## 1 - Simulate fledging and capture data (script 2)

In a study area, corresponding to the "sample area" of a capture site :

-   3-5 breeding pairs (capturable during all the study period)

-   Each couple lays at a time t (sample from a normal distribution N($\mu_{site}$, $\sigma_{site}$), $N_{eggs}$ sample from a Poisson distribution Pois($\lambda_{eggs}$).

-   Each egg has a probability $P_{fledge}$ to fledge which depends on the laying date.

-   Once fledge, survival of juveniles may evolve according to time [@naef-daenzer2001]

-   Still need to implement:

    -   Varying clutch size with laying date

    -   Decreasing fledging survival throught time

    -   Varying laying date variances between years

Simulating CES capture protocol would just correspond in sampling $N_{sessions}$ times (every 2 weeks for instance) all the available individuals (juveniles + adults).

Save data and apply model. (different version of the model)

## References

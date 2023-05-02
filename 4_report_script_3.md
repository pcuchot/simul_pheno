4_report_script_3
================
Paul Cuchot
02/05/2023

## Test model with simulated data

### Library

``` r
# bayesian modeling + plot
require(R2jags)
```

    ## Le chargement a nécessité le package : R2jags

    ## Warning: le package 'R2jags' a été compilé avec la version R 4.1.3

    ## Le chargement a nécessité le package : rjags

    ## Le chargement a nécessité le package : coda

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

    ## 
    ## Attachement du package : 'R2jags'

    ## L'objet suivant est masqué depuis 'package:coda':
    ## 
    ##     traceplot

``` r
require(MASS)
```

    ## Le chargement a nécessité le package : MASS

``` r
require(mcmcplots)
```

    ## Le chargement a nécessité le package : mcmcplots

    ## Warning: le package 'mcmcplots' a été compilé avec la version R 4.1.3

    ## Registered S3 method overwritten by 'mcmcplots':
    ##   method        from  
    ##   as.mcmc.rjags R2jags

``` r
# general
require(tidyverse)
```

    ## Le chargement a nécessité le package : tidyverse

    ## Warning: le package 'tidyverse' a été compilé avec la version R 4.1.3

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.4.0     v purrr   0.3.4
    ## v tibble  3.1.6     v dplyr   1.0.8
    ## v tidyr   1.2.0     v stringr 1.4.0
    ## v readr   2.1.2     v forcats 0.5.1

    ## Warning: le package 'ggplot2' a été compilé avec la version R 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()
    ## x dplyr::select() masks MASS::select()

``` r
require(bayesplot)
```

    ## Le chargement a nécessité le package : bayesplot

    ## This is bayesplot version 1.8.1

    ## - Online documentation and vignettes at mc-stan.org/bayesplot

    ## - bayesplot theme set to bayesplot::theme_default()

    ##    * Does _not_ affect other ggplot2 plots

    ##    * See ?bayesplot_theme_set for details on theme setting

### Simulate data

``` r
# load function 
source("2_function_sim_data.R")
# number of year to simlate
Nyears = 10

prod <- simul_data(
  # 5 pairs of breeders per year
  n_breeders = 5,
  # Nyears
  n_years = Nyears,
  # CES start (julian days)
  start_ces = 100,
  # CES end (julian days)
  end_ces = 220,
  # sessions per year 
  n_session = 5,
  # mean laying date for this site 
  mean_ld_site = 120)
```

### plot productivity through time (sim data)

![](4_report_script_3_files/figure-gfm/pressure-1.png)<!-- -->

### Structure data for JAGS

``` r
prod_f <- prod$capt_sess %>%
  mutate(an = as.numeric(year))

# data for the model
data <- list(nt = prod_f$n_capt_juveniles+prod_f$n_capt_adults,
             n0 = prod_f$n_capt_juveniles,
             date = as.numeric(prod_f$t),
             N = nrow(prod_f),
             N_an = length(unique(prod_f$an)),
             an = prod_f$an
)
```

### Define model

``` r
model <- "model{
 
  # loop on capture session
 
  for(i in 1:N){
    
    ## likelihood
    n0[i] ~ dbin(p[i], nt[i])
    
    p[i] <- asig[an[i]]/(1+exp((csig[an[i]]-date[i])/dsig[an[i]]))
  }
 
  # loop on an

  for (ii in 1:N_an){
    
    # asig[ii] <- alpha_asig + random_asig_an[an_n[ii]]
    # csig[ii] <- alpha_csig + random_csig_an[an_n[ii]]
    # dsig[ii] <- alpha_dsig + random_dsig_an[an_n[ii]]
    
    asig[ii] <- alpha_asig + random_asig_an[ii]
    csig[ii] <- alpha_csig + random_csig_an[ii]
    dsig[ii] <- alpha_dsig + random_dsig_an[ii]
    
    random_asig_an[ii] ~ dnorm(0, tau_asig_an)
    random_csig_an[ii] ~ dnorm(0, tau_csig_an)
    random_dsig_an[ii] ~ dnorm(0, tau_dsig_an)
 
  }

  alpha_asig ~ dnorm(0.7,0.1)T(0,1)
  alpha_csig ~ dnorm(150,0.1)T(100,200)
  alpha_dsig ~ dnorm(7,0.1)T(0,15)
 
#folded Cauchy version
  sigma_asig_an ~ dt(0, 0.01, 1)
  tau_asig_an <- pow(sigma_asig_an, -2)
 
#folded Cauchy version
  sigma_csig_an ~ dt(0, 0.01, 1)
  tau_csig_an <- pow(sigma_csig_an, -2)
  
  sigma_dsig_an ~ dt(0, 0.01, 1)
  tau_dsig_an <- pow(sigma_dsig_an, -2)

}
"
```

### Create initial values and parameters to save

``` r
init_f <-  function(){
  
  list(
    
    alpha_asig = runif(1, 0, 1),
    alpha_dsig = runif(1, 3, 10),
    alpha_csig = rnorm(1, 150, 20),
    
    sigma_asig_an = rnorm(1, 3, 20),
    sigma_csig_an = rnorm(1, 3, 20),
    sigma_dsig_an = rnorm(1, 3, 20),
    
    random_asig_an = rnorm(Nyears, 0, 0.001),
    random_csig_an = rnorm(Nyears, 0, 0.001),
    random_dsig_an = rnorm(Nyears, 0, 0.001)
    
  )
}

inits <- list(init1 = init_f(), init2=init_f(), init3=init_f())

# parameter to save 
parameters3 <- c("asig","csig","dsig",
                 "random_asig_an","random_csig_an","random_dsig_an",
                 "sigma_csig_an","sigma_dsig_an","sigma_asig_an")
```

### Run model

``` r
# run model 
md_1 <- jags(data = data,
             parameters.to.save = parameters3,
             model.file = textConnection(model),
             inits = inits,
             n.chains = 3,
             n.iter = 10000,
             n.burnin = 3000)
```

    ## module glm loaded

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 50
    ##    Unobserved stochastic nodes: 36
    ##    Total graph size: 533
    ## 
    ## Initializing model

### Look at convergence

``` r
md_1
```

    ## Inference for Bugs model at "4", fit using jags,
    ##  3 chains, each with 10000 iterations (first 3000 discarded), n.thin = 7
    ##  n.sims = 3000 iterations saved
    ##                    mu.vect sd.vect    2.5%     25%     50%     75%   97.5%
    ## asig[1]              0.797   0.042   0.705   0.773   0.799   0.824   0.874
    ## asig[2]              0.796   0.040   0.709   0.772   0.798   0.823   0.867
    ## asig[3]              0.807   0.038   0.730   0.783   0.808   0.830   0.880
    ## asig[4]              0.808   0.043   0.723   0.781   0.807   0.833   0.898
    ## asig[5]              0.783   0.045   0.679   0.759   0.789   0.813   0.852
    ## asig[6]              0.800   0.041   0.715   0.776   0.802   0.827   0.879
    ## asig[7]              0.804   0.043   0.721   0.778   0.804   0.829   0.890
    ## asig[8]              0.804   0.042   0.717   0.778   0.804   0.829   0.888
    ## asig[9]              0.806   0.039   0.724   0.782   0.806   0.830   0.885
    ## asig[10]             0.798   0.039   0.715   0.776   0.800   0.823   0.869
    ## csig[1]            171.698   6.899 161.211 166.310 170.659 176.550 186.157
    ## csig[2]            147.407   6.712 134.197 142.215 147.831 152.863 158.414
    ## csig[3]            119.807   5.772 107.024 115.939 120.717 124.372 128.354
    ## csig[4]            171.797   6.277 162.067 166.950 170.916 176.119 185.233
    ## csig[5]            126.755   5.324 111.856 125.050 128.145 129.884 134.315
    ## csig[6]            146.889   6.451 134.222 142.148 147.339 152.046 157.537
    ## csig[7]            150.936   7.567 134.819 145.548 152.465 157.062 161.610
    ## csig[8]            150.701   7.104 135.531 145.449 152.135 156.595 160.348
    ## csig[9]            146.224   6.027 134.788 141.786 146.418 150.732 157.063
    ## csig[10]           146.939   6.698 133.746 142.012 147.195 152.328 157.996
    ## dsig[1]              2.554   1.605   0.291   1.378   2.315   3.367   6.306
    ## dsig[2]              2.566   1.524   0.320   1.438   2.352   3.426   6.144
    ## dsig[3]              2.410   1.440   0.281   1.336   2.251   3.246   5.754
    ## dsig[4]              2.425   1.434   0.263   1.354   2.218   3.282   5.655
    ## dsig[5]              3.021   1.856   0.504   1.713   2.713   3.890   7.778
    ## dsig[6]              2.471   1.471   0.299   1.406   2.288   3.318   5.968
    ## dsig[7]              2.817   1.712   0.352   1.580   2.544   3.712   7.066
    ## dsig[8]              2.787   1.672   0.392   1.584   2.574   3.650   6.757
    ## dsig[9]              2.432   1.377   0.313   1.412   2.271   3.239   5.588
    ## dsig[10]             2.564   1.567   0.320   1.431   2.346   3.424   6.186
    ## random_asig_an[1]   -0.003   0.035  -0.083  -0.017  -0.001   0.012   0.067
    ## random_asig_an[2]   -0.004   0.033  -0.084  -0.016  -0.001   0.010   0.061
    ## random_asig_an[3]    0.006   0.031  -0.054  -0.008   0.002   0.020   0.079
    ## random_asig_an[4]    0.008   0.036  -0.060  -0.009   0.002   0.021   0.099
    ## random_asig_an[5]   -0.018   0.037  -0.114  -0.033  -0.007   0.002   0.038
    ## random_asig_an[6]    0.000   0.033  -0.074  -0.013   0.000   0.014   0.071
    ## random_asig_an[7]    0.003   0.036  -0.073  -0.010   0.001   0.017   0.089
    ## random_asig_an[8]    0.004   0.035  -0.070  -0.011   0.000   0.017   0.086
    ## random_asig_an[9]    0.005   0.032  -0.059  -0.009   0.002   0.019   0.080
    ## random_asig_an[10]  -0.002   0.032  -0.074  -0.015   0.000   0.012   0.064
    ## random_csig_an[1]   22.353   7.242  10.302  16.822  21.557  27.244  37.540
    ## random_csig_an[2]   -1.937   7.088 -16.049  -7.083  -1.746   3.310  10.442
    ## random_csig_an[3]  -29.538   6.379 -43.332 -33.725 -29.005 -24.803 -18.662
    ## random_csig_an[4]   22.453   6.645  11.301  17.455  21.846  26.922  36.280
    ## random_csig_an[5]  -22.589   5.870 -38.030 -25.343 -21.695 -18.661 -13.628
    ## random_csig_an[6]   -2.455   6.827 -15.938  -7.450  -2.035   2.723   9.324
    ## random_csig_an[7]    1.592   7.824 -14.659  -4.054   2.788   7.558  13.968
    ## random_csig_an[8]    1.356   7.365 -13.970  -3.861   2.331   6.978  13.372
    ## random_csig_an[9]   -3.120   6.502 -15.836  -7.650  -3.111   1.495   9.280
    ## random_csig_an[10]  -2.406   6.981 -15.935  -7.358  -2.218   2.739   9.949
    ## random_dsig_an[1]   -0.122   1.175  -2.678  -0.555  -0.052   0.319   2.389
    ## random_dsig_an[2]   -0.110   1.118  -2.769  -0.516  -0.038   0.319   2.221
    ## random_dsig_an[3]   -0.266   1.164  -3.041  -0.681  -0.082   0.254   1.884
    ## random_dsig_an[4]   -0.251   1.137  -3.047  -0.679  -0.086   0.257   1.876
    ## random_dsig_an[5]    0.345   1.262  -1.716  -0.204   0.096   0.710   3.753
    ## random_dsig_an[6]   -0.205   1.117  -2.899  -0.611  -0.071   0.265   1.943
    ## random_dsig_an[7]    0.140   1.175  -2.053  -0.319   0.032   0.493   2.982
    ## random_dsig_an[8]    0.111   1.132  -2.207  -0.343   0.020   0.509   2.692
    ## random_dsig_an[9]   -0.245   1.131  -3.138  -0.594  -0.069   0.236   1.802
    ## random_dsig_an[10]  -0.112   1.174  -2.676  -0.542  -0.042   0.322   2.282
    ## sigma_asig_an       -0.005   0.041  -0.088  -0.030  -0.008   0.020   0.077
    ## sigma_csig_an        6.003  17.520 -25.005 -14.330  14.346  18.415  28.218
    ## sigma_dsig_an       -0.092   1.393  -2.854  -0.880  -0.120   0.674   2.793
    ## deviance            87.386   4.508  80.543  84.197  86.564  89.919  98.193
    ##                     Rhat n.eff
    ## asig[1]            1.001  2600
    ## asig[2]            1.001  2300
    ## asig[3]            1.001  2900
    ## asig[4]            1.001  3000
    ## asig[5]            1.005   700
    ## asig[6]            1.001  2100
    ## asig[7]            1.001  3000
    ## asig[8]            1.001  3000
    ## asig[9]            1.003   910
    ## asig[10]           1.001  2600
    ## csig[1]            1.001  3000
    ## csig[2]            1.001  3000
    ## csig[3]            1.001  3000
    ## csig[4]            1.001  3000
    ## csig[5]            1.007   850
    ## csig[6]            1.001  3000
    ## csig[7]            1.001  3000
    ## csig[8]            1.002  1100
    ## csig[9]            1.001  2900
    ## csig[10]           1.001  3000
    ## dsig[1]            1.005   460
    ## dsig[2]            1.004   520
    ## dsig[3]            1.004   830
    ## dsig[4]            1.004   740
    ## dsig[5]            1.007   320
    ## dsig[6]            1.007   460
    ## dsig[7]            1.005   470
    ## dsig[8]            1.010   330
    ## dsig[9]            1.005   410
    ## dsig[10]           1.009   240
    ## random_asig_an[1]  1.006   720
    ## random_asig_an[2]  1.001  3000
    ## random_asig_an[3]  1.002  3000
    ## random_asig_an[4]  1.002  1700
    ## random_asig_an[5]  1.006  1000
    ## random_asig_an[6]  1.002  3000
    ## random_asig_an[7]  1.004  2800
    ## random_asig_an[8]  1.002  1900
    ## random_asig_an[9]  1.003  1200
    ## random_asig_an[10] 1.001  3000
    ## random_csig_an[1]  1.003  2400
    ## random_csig_an[2]  1.001  3000
    ## random_csig_an[3]  1.001  3000
    ## random_csig_an[4]  1.001  3000
    ## random_csig_an[5]  1.005   630
    ## random_csig_an[6]  1.001  3000
    ## random_csig_an[7]  1.001  2500
    ## random_csig_an[8]  1.003   880
    ## random_csig_an[9]  1.001  3000
    ## random_csig_an[10] 1.001  2300
    ## random_dsig_an[1]  1.010  3000
    ## random_dsig_an[2]  1.008   600
    ## random_dsig_an[3]  1.007   870
    ## random_dsig_an[4]  1.005  2100
    ## random_dsig_an[5]  1.005  2000
    ## random_dsig_an[6]  1.008  3000
    ## random_dsig_an[7]  1.011  3000
    ## random_dsig_an[8]  1.001  2000
    ## random_dsig_an[9]  1.004  2300
    ## random_dsig_an[10] 1.018  1600
    ## sigma_asig_an      1.048    52
    ## sigma_csig_an      6.311     3
    ## sigma_dsig_an      1.067    37
    ## deviance           1.002  1700
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = var(deviance)/2)
    ## pD = 10.2 and DIC = 97.5
    ## DIC is an estimate of expected predictive error (lower deviance is better).

### quickly compare mean xmid with ‘real pheno’

``` r
data.frame(estim = md_1$BUGSoutput$mean$c, 
           real = prod$mean_ld_year$mean_ld, 
           year = 1:Nyears)%>%
  ggplot(aes(x = estim, y = real, label = year))+
  geom_text(vjust = 1.5)+
  geom_point()+
  theme_bw()
```

![](4_report_script_3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

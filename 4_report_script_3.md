4_report_script_3
================
Paul Cuchot
02/05/2023

## Test model with simulated data

### Library

``` r
# bayesian modeling + plot
require(R2jags)
require(MASS)
require(mcmcplots)

# general
require(tidyverse)
require(bayesplot)
```

### Simulate data

#### Assumption for simulated data

-   Same number of breeders between years (5 pairs)
-   Mean number of eggs = 8 (sd = 5)
-   All chicks fledge and survive
-   Adults survive during the breeding season
-   Between years, sessions are on the same days
-   Birds do not migrate or immigrate during the breeding season
-   Same capture probability between juveniles and adults

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
  # mean laying date for this site 
  mean_ld_site = 120,
  
  # CES start (julian days)
  start_ces = 100,
  # CES end (julian days)
  end_ces = 220,
  # sessions per year 
  n_session = 5)
```

### plot productivity through time (sim data)

Productivity throught breeding period for a single site. The vertical
lines represent the simulated mean laying for each year.

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
    ## asig[1]              0.774   0.047   0.675   0.743   0.774   0.804   0.872
    ## asig[2]              0.788   0.051   0.693   0.756   0.784   0.819   0.897
    ## asig[3]              0.739   0.059   0.595   0.708   0.746   0.778   0.833
    ## asig[4]              0.776   0.048   0.684   0.747   0.776   0.805   0.874
    ## asig[5]              0.762   0.047   0.657   0.735   0.766   0.793   0.849
    ## asig[6]              0.761   0.051   0.651   0.732   0.765   0.794   0.854
    ## asig[7]              0.767   0.054   0.649   0.736   0.769   0.800   0.871
    ## asig[8]              0.745   0.055   0.613   0.714   0.752   0.782   0.838
    ## asig[9]              0.772   0.047   0.678   0.742   0.772   0.801   0.866
    ## asig[10]             0.763   0.049   0.659   0.733   0.765   0.795   0.855
    ## csig[1]            147.769   4.357 137.907 145.378 148.160 150.635 155.441
    ## csig[2]            147.978   4.344 138.200 145.543 148.305 150.849 155.590
    ## csig[3]            153.338   6.222 143.433 149.077 152.308 157.156 167.099
    ## csig[4]            147.840   4.395 137.612 145.326 148.277 150.800 155.463
    ## csig[5]            148.818   4.682 138.710 146.210 149.000 151.828 157.669
    ## csig[6]            143.396   7.047 128.529 138.508 145.181 148.815 153.557
    ## csig[7]            144.893   7.611 126.736 141.477 146.938 150.014 155.583
    ## csig[8]            151.236   5.216 141.393 147.848 150.839 154.404 161.962
    ## csig[9]            148.208   4.386 138.076 145.763 148.637 151.126 155.809
    ## csig[10]           149.131   4.643 139.257 146.324 149.240 152.154 158.044
    ## dsig[1]              4.013   2.520   0.332   2.128   3.746   5.475   9.797
    ## dsig[2]              4.365   2.699   0.255   2.325   4.051   5.934  10.444
    ## dsig[3]              8.274   4.173   1.667   5.231   7.773  10.685  17.637
    ## dsig[4]              4.315   2.823   0.288   2.195   3.921   5.842  11.226
    ## dsig[5]              5.085   3.021   0.489   2.898   4.675   6.772  11.897
    ## dsig[6]             10.397   4.630   3.095   7.137   9.910  13.015  21.223
    ## dsig[7]             12.183   4.699   3.166   9.375  12.154  14.984  21.957
    ## dsig[8]              7.297   3.680   1.290   4.804   6.810   9.310  16.039
    ## dsig[9]              4.705   2.772   0.369   2.718   4.384   6.298  11.003
    ## dsig[10]             5.638   3.332   0.600   3.259   5.177   7.436  13.270
    ## random_asig_an[1]    0.010   0.043  -0.074  -0.011   0.003   0.028   0.109
    ## random_asig_an[2]    0.024   0.047  -0.052  -0.003   0.012   0.045   0.141
    ## random_asig_an[3]   -0.025   0.050  -0.152  -0.048  -0.011   0.003   0.048
    ## random_asig_an[4]    0.012   0.043  -0.062  -0.008   0.004   0.030   0.118
    ## random_asig_an[5]   -0.002   0.040  -0.092  -0.019   0.000   0.016   0.084
    ## random_asig_an[6]   -0.003   0.044  -0.104  -0.021   0.000   0.017   0.084
    ## random_asig_an[7]    0.003   0.046  -0.097  -0.016   0.001   0.021   0.104
    ## random_asig_an[8]   -0.019   0.046  -0.136  -0.038  -0.008   0.005   0.061
    ## random_asig_an[9]    0.008   0.041  -0.072  -0.012   0.002   0.026   0.102
    ## random_asig_an[10]  -0.001   0.041  -0.090  -0.019   0.000   0.017   0.087
    ## random_csig_an[1]   -1.015   4.035 -10.529  -2.897  -0.423   1.010   6.634
    ## random_csig_an[2]   -0.806   4.138 -10.342  -2.824  -0.325   1.281   7.355
    ## random_csig_an[3]    4.554   6.041  -3.884   0.139   2.836   8.132  18.543
    ## random_csig_an[4]   -0.944   4.071 -10.703  -2.806  -0.357   1.170   7.035
    ## random_csig_an[5]    0.034   4.463  -9.823  -2.040   0.017   2.215   9.429
    ## random_csig_an[6]   -5.388   6.533 -20.027  -9.732  -3.320  -0.172   3.191
    ## random_csig_an[7]   -3.891   7.001 -21.700  -6.638  -1.437   0.276   5.704
    ## random_csig_an[8]    2.452   4.990  -6.572  -0.350   1.407   5.194  14.201
    ## random_csig_an[9]   -0.576   4.191 -10.083  -2.523  -0.123   1.543   7.864
    ## random_csig_an[10]   0.347   4.253  -8.733  -1.643   0.150   2.571   9.415
    ## random_dsig_an[1]   -2.521   2.800  -7.918  -4.427  -2.521  -0.561   2.913
    ## random_dsig_an[2]   -2.169   2.941  -7.871  -4.203  -2.044  -0.222   3.653
    ## random_dsig_an[3]    1.740   3.843  -4.748  -0.706   1.189   3.828  10.471
    ## random_dsig_an[4]   -2.220   3.012  -8.057  -4.230  -2.197  -0.207   3.868
    ## random_dsig_an[5]   -1.449   3.090  -7.150  -3.518  -1.446   0.328   5.034
    ## random_dsig_an[6]    3.862   4.109  -2.306   0.946   3.265   6.165  13.617
    ## random_dsig_an[7]    5.649   4.143  -0.894   2.848   5.313   8.043  14.429
    ## random_dsig_an[8]    0.763   3.512  -5.815  -1.371   0.474   2.566   8.692
    ## random_dsig_an[9]   -1.829   2.894  -7.638  -3.833  -1.723   0.076   4.040
    ## random_dsig_an[10]  -0.896   3.332  -6.943  -3.139  -0.881   0.975   6.157
    ## sigma_asig_an       -0.008   0.055  -0.116  -0.042  -0.010   0.028   0.101
    ## sigma_csig_an       -0.076   6.433 -12.421  -4.571   0.262   4.500  11.463
    ## sigma_dsig_an       -1.258   4.857  -9.117  -4.763  -2.821   3.188   7.829
    ## deviance           119.751   6.006 109.168 115.534 119.588 123.557 132.759
    ##                     Rhat n.eff
    ## asig[1]            1.004   560
    ## asig[2]            1.002  1100
    ## asig[3]            1.003   990
    ## asig[4]            1.001  3000
    ## asig[5]            1.001  3000
    ## asig[6]            1.005  3000
    ## asig[7]            1.002  3000
    ## asig[8]            1.003   990
    ## asig[9]            1.002  2200
    ## asig[10]           1.001  3000
    ## csig[1]            1.002  1100
    ## csig[2]            1.001  3000
    ## csig[3]            1.006   400
    ## csig[4]            1.001  3000
    ## csig[5]            1.001  3000
    ## csig[6]            1.003   680
    ## csig[7]            1.003   970
    ## csig[8]            1.003   880
    ## csig[9]            1.003   970
    ## csig[10]           1.002  3000
    ## dsig[1]            1.001  3000
    ## dsig[2]            1.002  2100
    ## dsig[3]            1.001  3000
    ## dsig[4]            1.001  3000
    ## dsig[5]            1.003  2100
    ## dsig[6]            1.002  1700
    ## dsig[7]            1.001  3000
    ## dsig[8]            1.003  3000
    ## dsig[9]            1.004  3000
    ## dsig[10]           1.001  3000
    ## random_asig_an[1]  1.007   310
    ## random_asig_an[2]  1.003   960
    ## random_asig_an[3]  1.005   850
    ## random_asig_an[4]  1.002  1600
    ## random_asig_an[5]  1.001  3000
    ## random_asig_an[6]  1.003  3000
    ## random_asig_an[7]  1.003  3000
    ## random_asig_an[8]  1.005  1100
    ## random_asig_an[9]  1.002  2500
    ## random_asig_an[10] 1.004  3000
    ## random_csig_an[1]  1.003   920
    ## random_csig_an[2]  1.004  3000
    ## random_csig_an[3]  1.007   360
    ## random_csig_an[4]  1.001  2800
    ## random_csig_an[5]  1.001  3000
    ## random_csig_an[6]  1.004   600
    ## random_csig_an[7]  1.003   840
    ## random_csig_an[8]  1.003   740
    ## random_csig_an[9]  1.004  1100
    ## random_csig_an[10] 1.005  2800
    ## random_dsig_an[1]  1.002  1600
    ## random_dsig_an[2]  1.002  1500
    ## random_dsig_an[3]  1.001  3000
    ## random_dsig_an[4]  1.001  2300
    ## random_dsig_an[5]  1.001  2500
    ## random_dsig_an[6]  1.003   820
    ## random_dsig_an[7]  1.001  3000
    ## random_dsig_an[8]  1.001  3000
    ## random_dsig_an[9]  1.001  3000
    ## random_dsig_an[10] 1.001  3000
    ## sigma_asig_an      1.161    17
    ## sigma_csig_an      1.023   100
    ## sigma_dsig_an      1.089    28
    ## deviance           1.003   900
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = var(deviance)/2)
    ## pD = 18.0 and DIC = 137.8
    ## DIC is an estimate of expected predictive error (lower deviance is better).

### quickly compare mean xmid with ‘real pheno’

``` r
data.frame(estim = md_1$BUGSoutput$mean$c, 
           real = prod$mean_ld_year$mean_ld, 
           year = 1:Nyears)%>%
  ggplot(aes(x = estim, y = real, label = year))+
  geom_text(vjust = 1.5)+
  geom_point()+
  theme_bw()+
  stat_smooth(method = lm)
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: The following aesthetics were dropped during statistical transformation: label
    ## i This can happen when ggplot fails to infer the correct grouping structure in
    ##   the data.
    ## i Did you forget to specify a `group` aesthetic or to convert a numerical
    ##   variable into a factor?

![](4_report_script_3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

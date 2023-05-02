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

The function used to simulate data is described at the end of this
document.

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
    ## asig[1]              0.710   0.071   0.551   0.665   0.717   0.761   0.829
    ## asig[2]              0.750   0.055   0.634   0.717   0.752   0.786   0.851
    ## asig[3]              0.659   0.086   0.475   0.603   0.668   0.725   0.795
    ## asig[4]              0.794   0.058   0.678   0.754   0.793   0.835   0.904
    ## asig[5]              0.802   0.061   0.689   0.760   0.798   0.841   0.926
    ## asig[6]              0.791   0.058   0.681   0.751   0.790   0.830   0.907
    ## asig[7]              0.767   0.058   0.646   0.731   0.768   0.805   0.879
    ## asig[8]              0.739   0.066   0.588   0.700   0.744   0.782   0.855
    ## asig[9]              0.789   0.060   0.675   0.750   0.787   0.828   0.911
    ## asig[10]             0.762   0.059   0.635   0.726   0.763   0.800   0.876
    ## csig[1]            155.055   7.259 138.196 150.820 156.799 160.066 166.294
    ## csig[2]            125.916   5.872 112.759 122.989 126.410 129.060 138.319
    ## csig[3]            135.438   7.612 121.541 130.929 134.261 139.319 153.520
    ## csig[4]            139.625   5.362 131.042 135.570 139.149 143.281 151.101
    ## csig[5]            147.946   6.287 135.476 143.423 148.375 152.854 158.637
    ## csig[6]            146.415   5.743 134.951 142.296 146.589 150.768 156.762
    ## csig[7]            145.907   5.893 134.580 141.624 146.109 150.259 156.541
    ## csig[8]            170.498   6.156 160.785 165.797 169.837 174.515 183.681
    ## csig[9]            146.243   6.010 134.465 141.965 146.465 150.793 156.894
    ## csig[10]           150.073   6.517 137.062 145.381 150.683 155.161 160.693
    ## dsig[1]              4.478   2.482   1.023   2.878   4.089   5.491  10.117
    ## dsig[2]              4.456   2.798   0.833   2.731   3.901   5.389  12.182
    ## dsig[3]              5.357   3.168   1.478   3.247   4.585   6.624  13.746
    ## dsig[4]              4.458   1.979   1.540   3.114   4.162   5.396   9.267
    ## dsig[5]              3.744   1.952   0.469   2.403   3.561   4.800   8.187
    ## dsig[6]              3.415   1.747   0.430   2.234   3.289   4.448   7.341
    ## dsig[7]              3.450   1.824   0.438   2.222   3.309   4.480   7.521
    ## dsig[8]              3.431   1.930   0.309   2.084   3.250   4.526   7.883
    ## dsig[9]              3.521   1.818   0.473   2.255   3.382   4.557   7.549
    ## dsig[10]             3.966   2.091   0.662   2.531   3.718   5.032   9.044
    ## random_asig_an[1]   -0.048   0.068  -0.201  -0.088  -0.037   0.000   0.066
    ## random_asig_an[2]   -0.008   0.055  -0.121  -0.039  -0.005   0.023   0.104
    ## random_asig_an[3]   -0.099   0.082  -0.289  -0.147  -0.085  -0.035   0.014
    ## random_asig_an[4]    0.037   0.060  -0.074  -0.001   0.030   0.075   0.166
    ## random_asig_an[5]    0.045   0.063  -0.065   0.001   0.036   0.084   0.186
    ## random_asig_an[6]    0.034   0.060  -0.076  -0.003   0.028   0.071   0.163
    ## random_asig_an[7]    0.010   0.060  -0.115  -0.023   0.009   0.045   0.133
    ## random_asig_an[8]   -0.019   0.063  -0.160  -0.052  -0.013   0.018   0.102
    ## random_asig_an[9]    0.032   0.063  -0.085  -0.006   0.025   0.069   0.168
    ## random_asig_an[10]   0.005   0.059  -0.119  -0.026   0.003   0.038   0.127
    ## random_csig_an[1]    6.363   7.405 -10.746   2.089   7.642  11.399  18.320
    ## random_csig_an[2]  -22.776   6.279 -36.420 -26.345 -22.520 -18.857 -10.762
    ## random_csig_an[3]  -13.254   7.782 -27.794 -17.927 -14.003  -9.126   4.279
    ## random_csig_an[4]   -9.067   5.751 -19.285 -13.149  -9.401  -5.402   3.057
    ## random_csig_an[5]   -0.745   6.549 -13.392  -5.320  -0.443   4.104  10.829
    ## random_csig_an[6]   -2.277   6.057 -14.274  -6.379  -2.085   2.173   8.864
    ## random_csig_an[7]   -2.785   6.211 -14.940  -7.038  -2.616   1.744   8.405
    ## random_csig_an[8]   21.806   6.344  11.082  17.196  21.334  25.994  34.775
    ## random_csig_an[9]   -2.449   6.320 -14.770  -7.001  -2.291   2.173   9.359
    ## random_csig_an[10]   1.381   6.639 -12.174  -3.224   1.798   6.396  13.212
    ## random_dsig_an[1]    0.355   2.005  -3.234  -0.489   0.093   1.032   4.797
    ## random_dsig_an[2]    0.334   2.126  -3.250  -0.619   0.042   0.982   5.883
    ## random_dsig_an[3]    1.234   2.424  -1.881  -0.080   0.511   1.996   7.517
    ## random_dsig_an[4]    0.335   1.635  -2.941  -0.373   0.135   1.017   3.977
    ## random_dsig_an[5]   -0.379   1.768  -4.746  -1.067  -0.129   0.434   3.041
    ## random_dsig_an[6]   -0.707   1.779  -5.339  -1.440  -0.308   0.200   2.354
    ## random_dsig_an[7]   -0.673   1.771  -4.761  -1.448  -0.306   0.212   2.409
    ## random_dsig_an[8]   -0.692   1.834  -5.222  -1.561  -0.306   0.244   2.703
    ## random_dsig_an[9]   -0.602   1.778  -4.792  -1.409  -0.244   0.287   2.603
    ## random_dsig_an[10]  -0.157   1.778  -4.080  -0.900  -0.034   0.567   3.811
    ## sigma_asig_an        0.008   0.091  -0.162  -0.068   0.024   0.078   0.161
    ## sigma_csig_an       13.939   4.010   7.935  11.151  13.346  16.078  23.570
    ## sigma_dsig_an        0.104   2.346  -4.419  -1.335   0.211   1.521   4.630
    ## deviance           108.070   6.709  96.385 103.281 107.588 112.354 122.075
    ##                     Rhat n.eff
    ## asig[1]            1.002  1400
    ## asig[2]            1.002  1700
    ## asig[3]            1.006   400
    ## asig[4]            1.002  1600
    ## asig[5]            1.003   860
    ## asig[6]            1.003   700
    ## asig[7]            1.001  3000
    ## asig[8]            1.001  3000
    ## asig[9]            1.006   370
    ## asig[10]           1.002  1900
    ## csig[1]            1.001  2100
    ## csig[2]            1.002  1700
    ## csig[3]            1.005   490
    ## csig[4]            1.001  3000
    ## csig[5]            1.002  2000
    ## csig[6]            1.001  3000
    ## csig[7]            1.001  3000
    ## csig[8]            1.002  1400
    ## csig[9]            1.001  3000
    ## csig[10]           1.001  3000
    ## dsig[1]            1.001  3000
    ## dsig[2]            1.005  1800
    ## dsig[3]            1.005   840
    ## dsig[4]            1.003   790
    ## dsig[5]            1.002  3000
    ## dsig[6]            1.004  2200
    ## dsig[7]            1.006  3000
    ## dsig[8]            1.002  3000
    ## dsig[9]            1.009  1900
    ## dsig[10]           1.002  1800
    ## random_asig_an[1]  1.003   890
    ## random_asig_an[2]  1.003  1100
    ## random_asig_an[3]  1.007   320
    ## random_asig_an[4]  1.001  2000
    ## random_asig_an[5]  1.002  1300
    ## random_asig_an[6]  1.002  1000
    ## random_asig_an[7]  1.001  3000
    ## random_asig_an[8]  1.001  3000
    ## random_asig_an[9]  1.004   510
    ## random_asig_an[10] 1.001  3000
    ## random_csig_an[1]  1.002  1500
    ## random_csig_an[2]  1.002  1200
    ## random_csig_an[3]  1.005   430
    ## random_csig_an[4]  1.001  3000
    ## random_csig_an[5]  1.002  1300
    ## random_csig_an[6]  1.001  2900
    ## random_csig_an[7]  1.001  3000
    ## random_csig_an[8]  1.002  2000
    ## random_csig_an[9]  1.001  3000
    ## random_csig_an[10] 1.001  2100
    ## random_dsig_an[1]  1.004  1700
    ## random_dsig_an[2]  1.008  1500
    ## random_dsig_an[3]  1.005   940
    ## random_dsig_an[4]  1.002  1300
    ## random_dsig_an[5]  1.004  2800
    ## random_dsig_an[6]  1.008   480
    ## random_dsig_an[7]  1.009   450
    ## random_dsig_an[8]  1.005  1500
    ## random_dsig_an[9]  1.010   410
    ## random_dsig_an[10] 1.006  3000
    ## sigma_asig_an      1.144    19
    ## sigma_csig_an      1.002  1400
    ## sigma_dsig_an      1.058    42
    ## deviance           1.005   440
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = var(deviance)/2)
    ## pD = 22.4 and DIC = 130.5
    ## DIC is an estimate of expected predictive error (lower deviance is better).

### Quickly compare mean xmid with ‘real pheno’

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

## Annexe

### Function to simulate data

``` r
simul_data <- function(n_breeders,
                       n_years, 
                       n_session, 
                       start_ces = 80,
                       end_ces = 120,
                       mean_ld_site = 120){
  
  # mean laying date (among breeding individuals - change between years)
  mean_ld <- round(rnorm(n_years, mean_ld_site, 15)) # (real pheno)
  
  sd_ld <- 7
  
  # number of eggs per pair
  mean_eggs <- 8
  sd_eggs <- 5
  
  # final data_set
  df_site <- data.frame(t = NA,
                        n_capt_adults = NA,
                        n_capt_juveniles = NA, 
                        prod = NA,
                        year = NA)[0,] 
  
  for(k in 1:n_years){
    
    # sample n_breeders laying events
    ld_dates <- round(rnorm(n_breeders, 
                            mean = mean_ld[k], 
                            sd = sd_ld))
    
    # fledglings dates 
    fledgl_dates <- round(ld_dates)+35
    
    # number of eggs per pair
    n_eggs <- abs(round(rnorm(n_breeders, mean_eggs, sd_eggs)))
    
    # create a dataframe (one row per breeding pair)
    df_breed <- data.frame(
      ld_dates = ld_dates,
      n_eggs = n_eggs,
      fledgl_dates = fledgl_dates
    )
    
    # sample (as CES design) 
    
    # choose days for capture session
    t_capt <- round(seq(start_ces, end_ces, 
                        length.out = n_session))
    
    # Number of sample individual per session (~capture effort)
    mean_n_capt <- 9
    
    # Dataframe with n_adults and n_juveniles captured per session
    df_session <- data.frame(t = t_capt,
                             n_capt_adults = NA,
                             n_capt_juveniles = NA,
                             prod = NA,
                             year = as.character(k)) 
    
    
    for(i in t_capt){
      
      # catchable adults
      n_adults <- n_breeders * 2
      
      # catchable juveniles
      n_juveniles <- sum(df_breed[df_breed$fledgl_dates < i ,]$n_eggs)
      
      # sample birds among available individuals
      capt_indiv <- sample( # how many individuals captured
        c(rep(0,n_adults), 
          rep(1,n_juveniles)),
        round(rnorm(1,mean_n_capt)), replace = TRUE)
      
      # how many adults
      df_session[df_session$t ==i,"n_capt_adults"] <- sum(capt_indiv == 0)
      # how many juveniles
      df_session[df_session$t ==i,"n_capt_juveniles"] <- sum(capt_indiv == 1) 
      
    }
    
    # calculate productivity
    df_session <- df_session%>%
      mutate(prod = n_capt_juveniles/(n_capt_adults+n_capt_juveniles))
    
    df_site <- rbind(df_site, df_session)
  }
  
  # df with mean ld per year ("real breeding date")
  df_mean_ld <- data.frame(year = 1:n_years,
                           mean_ld = mean_ld)
  
  return(list(capt_sess = df_site, 
              mean_ld_year = df_mean_ld))
  
}
```

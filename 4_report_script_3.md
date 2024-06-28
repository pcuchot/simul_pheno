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

- Same number of breeders between years (5 pairs)
- Mean number of eggs = 8 (sd = 5)
- All chicks fledge and survive
- Adults survive during the breeding season
- Between years, sessions are on the same days
- Birds do not migrate or immigrate during the breeding season
- Same capture probability between juveniles and adults

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
  # 10 years
  n_years = 10,
  # CES start (julian days)
  start_ces = 120,
  # CES end (julian days)
  end_ces = 200,
  # sessions per year
  n_session = 9,
  # mean laying date for this site
  mean_ld_site = 120)
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
    ##    Observed stochastic nodes: 90
    ##    Unobserved stochastic nodes: 36
    ##    Total graph size: 893
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
    ## asig[1]              0.751   0.043   0.666   0.725   0.750   0.777   0.840
    ## asig[2]              0.739   0.047   0.640   0.714   0.740   0.768   0.829
    ## asig[3]              0.724   0.051   0.598   0.697   0.729   0.757   0.807
    ## asig[4]              0.730   0.056   0.603   0.704   0.736   0.764   0.826
    ## asig[5]              0.748   0.040   0.666   0.723   0.748   0.774   0.829
    ## asig[6]              0.759   0.040   0.685   0.731   0.756   0.783   0.847
    ## asig[7]              0.735   0.040   0.651   0.712   0.737   0.761   0.809
    ## asig[8]              0.727   0.045   0.624   0.701   0.731   0.757   0.808
    ## asig[9]              0.755   0.042   0.676   0.727   0.752   0.779   0.846
    ## asig[10]             0.750   0.039   0.674   0.724   0.749   0.775   0.830
    ## csig[1]            174.290   2.174 170.273 172.719 174.255 175.870 178.410
    ## csig[2]            176.895   2.807 171.274 174.958 177.154 178.949 181.938
    ## csig[3]            178.027   3.259 172.051 175.865 178.097 179.999 184.874
    ## csig[4]            182.365   3.486 175.686 180.178 181.713 184.591 189.704
    ## csig[5]            149.399   2.457 143.437 148.196 149.650 150.883 153.811
    ## csig[6]            147.048   2.420 142.001 145.413 147.312 148.843 151.244
    ## csig[7]            156.982   2.489 151.768 155.353 157.303 158.764 161.199
    ## csig[8]            139.685   3.959 132.555 137.083 139.286 141.900 148.435
    ## csig[9]            174.593   2.144 170.689 172.957 174.619 176.131 178.567
    ## csig[10]           138.615   2.668 132.263 137.266 139.025 140.178 143.317
    ## dsig[1]              1.406   0.903   0.105   0.752   1.270   1.889   3.626
    ## dsig[2]              2.077   1.254   0.230   1.173   1.843   2.757   5.112
    ## dsig[3]              2.217   1.394   0.247   1.231   1.954   2.929   5.695
    ## dsig[4]              2.936   1.611   0.540   1.630   2.744   4.036   6.409
    ## dsig[5]              2.254   1.186   0.480   1.379   2.066   2.888   5.030
    ## dsig[6]              1.899   1.126   0.230   1.107   1.698   2.512   4.602
    ## dsig[7]              1.887   1.219   0.209   1.028   1.650   2.479   4.767
    ## dsig[8]              2.879   1.888   0.380   1.401   2.373   4.047   7.277
    ## dsig[9]              1.392   0.871   0.111   0.744   1.275   1.878   3.427
    ## dsig[10]             2.451   1.400   0.472   1.421   2.165   3.182   5.848
    ## random_asig_an[1]    0.010   0.038  -0.062  -0.008   0.004   0.025   0.099
    ## random_asig_an[2]   -0.003   0.039  -0.095  -0.017  -0.001   0.014   0.078
    ## random_asig_an[3]   -0.018   0.041  -0.120  -0.034  -0.007   0.004   0.050
    ## random_asig_an[4]   -0.011   0.043  -0.116  -0.026  -0.003   0.009   0.064
    ## random_asig_an[5]    0.006   0.036  -0.065  -0.011   0.002   0.021   0.091
    ## random_asig_an[6]    0.017   0.038  -0.045  -0.004   0.008   0.035   0.112
    ## random_asig_an[7]   -0.006   0.033  -0.083  -0.020  -0.003   0.010   0.061
    ## random_asig_an[8]   -0.015   0.037  -0.103  -0.032  -0.006   0.005   0.049
    ## random_asig_an[9]    0.013   0.038  -0.053  -0.006   0.005   0.029   0.107
    ## random_asig_an[10]   0.008   0.035  -0.056  -0.009   0.003   0.023   0.092
    ## random_csig_an[1]   21.800   3.611  15.115  19.351  21.655  24.216  29.057
    ## random_csig_an[2]   24.404   3.988  16.806  21.677  24.361  27.045  32.418
    ## random_csig_an[3]   25.536   4.409  17.110  22.559  25.435  28.404  34.714
    ## random_csig_an[4]   29.875   4.528  21.431  26.861  29.518  32.801  39.333
    ## random_csig_an[5]   -3.092   3.719 -10.565  -5.602  -3.046  -0.668   4.171
    ## random_csig_an[6]   -5.442   3.686 -12.516  -7.958  -5.459  -2.933   1.757
    ## random_csig_an[7]    4.492   3.745  -2.764   1.841   4.500   6.986  12.155
    ## random_csig_an[8]  -12.806   4.812 -21.223 -16.097 -13.156  -9.908  -2.238
    ## random_csig_an[9]   22.103   3.565  15.158  19.758  22.046  24.516  29.345
    ## random_csig_an[10] -13.876   3.777 -21.614 -16.304 -13.820 -11.481  -6.446
    ## random_dsig_an[1]   -0.830   1.121  -3.465  -1.504  -0.566  -0.017   0.792
    ## random_dsig_an[2]   -0.160   1.061  -2.570  -0.659  -0.060   0.348   2.083
    ## random_dsig_an[3]   -0.019   1.099  -2.303  -0.536  -0.017   0.463   2.432
    ## random_dsig_an[4]    0.700   1.104  -1.033  -0.011   0.432   1.285   3.374
    ## random_dsig_an[5]    0.018   0.988  -2.155  -0.426   0.005   0.456   2.159
    ## random_dsig_an[6]   -0.338   1.032  -2.709  -0.902  -0.135   0.188   1.619
    ## random_dsig_an[7]   -0.349   1.089  -2.784  -0.871  -0.168   0.165   1.763
    ## random_dsig_an[8]    0.642   1.358  -1.445  -0.136   0.247   1.285   4.074
    ## random_dsig_an[9]   -0.844   1.122  -3.428  -1.501  -0.554  -0.021   0.793
    ## random_dsig_an[10]   0.215   1.087  -1.899  -0.286   0.068   0.669   2.778
    ## sigma_asig_an        0.008   0.046  -0.086  -0.021   0.010   0.037   0.097
    ## sigma_csig_an        6.639  19.549 -27.956 -16.499  16.322  20.546  30.496
    ## sigma_dsig_an        0.249   1.438  -2.709  -0.670   0.362   1.218   2.752
    ## deviance           158.746   6.203 148.424 154.307 158.156 162.520 172.207
    ##                     Rhat n.eff
    ## asig[1]            1.002  2700
    ## asig[2]            1.006   740
    ## asig[3]            1.005   450
    ## asig[4]            1.008   540
    ## asig[5]            1.001  3000
    ## asig[6]            1.002  1600
    ## asig[7]            1.004   600
    ## asig[8]            1.006   430
    ## asig[9]            1.004  1200
    ## asig[10]           1.001  3000
    ## csig[1]            1.001  3000
    ## csig[2]            1.001  3000
    ## csig[3]            1.001  3000
    ## csig[4]            1.001  3000
    ## csig[5]            1.001  3000
    ## csig[6]            1.001  3000
    ## csig[7]            1.001  3000
    ## csig[8]            1.004   600
    ## csig[9]            1.001  3000
    ## csig[10]           1.001  3000
    ## dsig[1]            1.001  3000
    ## dsig[2]            1.002  3000
    ## dsig[3]            1.001  3000
    ## dsig[4]            1.007   500
    ## dsig[5]            1.001  3000
    ## dsig[6]            1.001  3000
    ## dsig[7]            1.001  3000
    ## dsig[8]            1.002  1300
    ## dsig[9]            1.003  1400
    ## dsig[10]           1.005   760
    ## random_asig_an[1]  1.004  2500
    ## random_asig_an[2]  1.005  2300
    ## random_asig_an[3]  1.005   580
    ## random_asig_an[4]  1.008  1200
    ## random_asig_an[5]  1.003  3000
    ## random_asig_an[6]  1.009   450
    ## random_asig_an[7]  1.006  2700
    ## random_asig_an[8]  1.006  1000
    ## random_asig_an[9]  1.014   380
    ## random_asig_an[10] 1.005  2100
    ## random_csig_an[1]  1.003   830
    ## random_csig_an[2]  1.001  3000
    ## random_csig_an[3]  1.001  2100
    ## random_csig_an[4]  1.002  1400
    ## random_csig_an[5]  1.002  1500
    ## random_csig_an[6]  1.003   860
    ## random_csig_an[7]  1.003   870
    ## random_csig_an[8]  1.004   640
    ## random_csig_an[9]  1.002  1300
    ## random_csig_an[10] 1.002  1300
    ## random_dsig_an[1]  1.002  1000
    ## random_dsig_an[2]  1.002  1500
    ## random_dsig_an[3]  1.001  3000
    ## random_dsig_an[4]  1.002  1100
    ## random_dsig_an[5]  1.001  2600
    ## random_dsig_an[6]  1.001  3000
    ## random_dsig_an[7]  1.001  3000
    ## random_dsig_an[8]  1.001  2000
    ## random_dsig_an[9]  1.002  1200
    ## random_dsig_an[10] 1.001  3000
    ## sigma_asig_an      1.018   160
    ## sigma_csig_an      6.744     3
    ## sigma_dsig_an      1.055    41
    ## deviance           1.002  1700
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = var(deviance)/2)
    ## pD = 19.2 and DIC = 178.0
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

    ## Warning: The following aesthetics were dropped during statistical transformation: label.
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in
    ##   the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
    ##   variable into a factor?

![](4_report_script_3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Annexe

### Function to simulate data

``` r
simul_data 
```

    ## function (n_breeders, n_years, n_session, start_ces = 80, end_ces = 120, 
    ##     mean_ld_site = 90) 
    ## {
    ##     mean_ld <- round(rnorm(n_years, mean_ld_site, 15))
    ##     sd_ld <- 7
    ##     mean_eggs <- 8
    ##     sd_eggs <- 2
    ##     df_site <- data.frame(t = NA, n_capt_adults = NA, n_capt_juveniles = NA, 
    ##         prod = NA, year = NA)[0, ]
    ##     for (k in 1:n_years) {
    ##         ld_dates <- round(rnorm(n_breeders, mean = mean_ld[k], 
    ##             sd = sd_ld))
    ##         fledgl_dates <- round(ld_dates) + 50
    ##         n_eggs <- abs(round(rnorm(n_breeders, mean_eggs, sd_eggs)))
    ##         df_breed <- data.frame(ld_date = ld_dates, n_egg = n_eggs, 
    ##             fledgl_dates = fledgl_dates)
    ##         df_breed <- df_breed %>% mutate(rap = ld_date/(mean_ld[k] * 
    ##             0.9), n_eggs = round(n_eggs * (1/rap)), egg_lost = n_egg - 
    ##             n_eggs) %>% arrange(ld_date)
    ##         t_capt <- round(seq(start_ces, end_ces, length.out = n_session))
    ##         mean_n_capt <- 9
    ##         df_session <- data.frame(t = t_capt, n_capt_adults = NA, 
    ##             n_capt_juveniles = NA, prod = NA, year = as.character(k))
    ##         for (i in t_capt) {
    ##             n_adults <- n_breeders * 2
    ##             n_juveniles <- sum(df_breed[df_breed$fledgl_date < 
    ##                 i, ]$n_eggs)
    ##             capt_indiv <- sample(c(rep(0, n_adults), rep(1, n_juveniles)), 
    ##                 round(rnorm(1, mean_n_capt, sd = 3)), replace = TRUE)
    ##             df_session[df_session$t == i, "n_capt_adults"] <- sum(capt_indiv == 
    ##                 0)
    ##             df_session[df_session$t == i, "n_capt_juveniles"] <- sum(capt_indiv == 
    ##                 1)
    ##         }
    ##         df_session <- df_session %>% mutate(prod = n_capt_juveniles/(n_capt_adults + 
    ##             n_capt_juveniles))
    ##         df_site <- rbind(df_site, df_session)
    ##     }
    ##     df_mean_ld <- data.frame(year = 1:n_years, mean_ld = mean_ld)
    ##     return(list(capt_sess = df_site, mean_ld_year = df_mean_ld))
    ## }
    ## <bytecode: 0x000001b3c9680a20>


# -------------------------------------------------------------------------
# Title : 3_model_pheno
# Author : Paul Cuchot
# Date : 27/04/23
# note :
# -------------------------------------------------------------------------

# Library -----------------------------------------------------------------

# bayesian modeling + plot
require(R2jags)
require(MASS)
require(mcmcplots)

# general
require(tidyverse)
require(bayesplot)

# simulate data
source("2_function_sim_data.R")

Nyears = 15

prod <- simul_data(
  # 5 pairs of breeders per year
  n_breeders = 5,
  # Nyears
  n_years = Nyears,
  # CES start (julian days)
  start_ces = 80,
  # CES end (julian days)
  end_ces = 220,
  # sessions per year 
  n_session = 8,
  # mean laying date for this site 
  mean_ld_site = 120)

# plot simulated
prod$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add laying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld, 
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  theme_light()

# Structure data ----------------------------------------------------------

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

# model3 ------------------------------------------------------------------

sink("model_simul3")
cat("
    model{
 
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
", fill=TRUE)
sink()

# print()


# Initial values ----------------------------------------------------------


init_f <-  function(){
  
  list(
    
    # asig = runif(Nyears, 0.5, 1),
    # dsig = runif(Nyears, 3,10),
    # csig = dnorm(Nyears, 150, 20),
    
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



# Run models --------------------------------------------------------------

# parameter to save md1
# parameters <- c("asig","csig","dsig",
#                 "sigma_csig","sigma_dsig","sigma_asig",
#                 "mu_csig","mu_dsig","mu_asig",
#                 "deviance"
# )

# parameter to save md2
# parameters2 <- c("asig","csig","dsig",
#                 "mu_csig","mu_dsig","mu_asig", 
#                 "deviance")

# parameter to save md3
parameters3 <- c("asig","csig","dsig",
                 "random_asig_an","random_csig_an","random_dsig_an",
                 "sigma_csig_an","sigma_dsig_an","sigma_asig_an")

# run model 
md_1 <- jags(data = data,
             parameters.to.save = parameters3,
             model.file = "model_simul3",
             inits = inits,
             n.chains = 3,
             n.iter = 10000,
             n.burnin = 3000)

# mcmcplot(md_1)


# quickly compare mean xmid with 'real pheno'
df_compare <- data.frame(estim = md_1$BUGSoutput$mean$c, 
           real = prod$mean_ld_year$mean_ld, 
           year = 1:Nyears)

# plot correlation between real and estimated phenology
df_compare%>%
  ggplot(aes(x = estim, y = real, label = year))+
  geom_text(vjust = 1.5)+
  geom_point()+
  theme_bw()







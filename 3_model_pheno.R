
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

# simulate data
source("2_function_sim_data.R")

prod <- simul_data(
  # 5 pairs of breeders per year
  n_breeders = 5,
  # 10 years
  n_years = 10,
  # CES start (julian days)
  start_ces = 80,
  # CES end (julian days)
  end_ces = 200,
  # sessions per year 
  n_session = 9,
  # mean laying date for this site 
  mean_ld_site = 120)

prod_f <- 
  mutate()

# Structure data ----------------------------------------------------------

prod_f <- prod$capt_sess %>%
  mutate(sitename = "1",
         species = "sp1",
         an_f = as.factor(year),
         an_n = as.numeric(an_f),
         site_f = as.character(sitename),
         site_n = as.numeric(site_f),
         sp_f = as.factor(species),
         sp_n = as.numeric(sp_f),
         an_sp = as.factor(paste(an_f,sp_f)),
         an_sp_n = as.numeric(an_sp)
         )

df_rd <- prod_f%>%
  distinct(an_sp_n, .keep_all = TRUE)

# data for the model
data <- list(nt = prod_f$n_capt_juveniles+prod_f$n_capt_adults,
             n0 = prod_f$n_capt_juveniles,
             date = as.numeric(prod_f$t),
             N = nrow(prod_f),
             N_an = length(unique(prod_f$an_n)),
             N_site = length(unique(prod_f$site_n)),
             an = prod_f$an_n,
             sp = as.numeric(as.factor(as.character(prod_f$sp_n))),
             an2 = df_rd$an_n,
             sp2 = as.numeric(as.factor(as.character(df_rd$sp_n))),
             site = df_rd$site_n,
             N_sp = n_distinct(df_rd$species),
             N_an_sp_n = n_distinct(prod_f$an_sp_n),
             an_sp_n = prod_f$an_sp_n,
             rd_an_sp_n = df_rd$an_sp_n
)

# Model 1 -----------------------------------------------------------------

sink("model_simul")
cat("
    model{

  # loop on capture session

  for(i in 1:N){ # session data

    ## likelihood
    n0[i] ~ dbin(p[i], nt[i])
    
    p[i] <- asig[an_sp_n[i]]/(1+exp((csig[an_sp_n[i]]-date[i])/dsig[an_sp_n[i]]))
  
  }


  # loop on year and species

  for (ii in 1:N_an_sp_n){
    

  # csig parameter
  
    csig[ii] ~ dnorm(mu[ii], tau_res_csig[sp2[ii]]) 

    mu[ii] <- c[ii] + random_csig_site[site[ii]]  

    
  # scale parameter
  
    dsig[ii] ~ dnorm(mu_dsig[ii], tau_res_dsig[sp2[ii]])
    
    mu_dsig[ii] <- d[ii] + random_dsig_site[site[ii]]


  # asymptote parameter  
  
    asig[ii] ~ dnorm(mu_asig[ii], tau_res_asig[sp2[ii]])T(0.3,1)
    
    mu_asig[ii] <- a[ii] + random_asig_site[site[ii]]
    
     
  c[ii] ~ dnorm(150,0.01)

  a[ii] ~ dnorm(0,0.01)T(0,1)

  d[ii] ~ dnorm(0,0.01)T(0,10)
    
    
  }
  

  # random effect site
  
  for(z in 1:N_site){ #number of sites
 
  # random site 
    random_csig_site[z] ~ dnorm(0, tau_csig_site)
    random_asig_site[z] ~ dnorm(0, tau_asig_site)
    random_dsig_site[z] ~ dnorm(0, tau_dsig_site)
    
  }
  
  for(s in 1:N_sp){ # for  variance estimation (per species)
  
    sigma_res_csig[s] ~ dt(0, 0.01, 1)T(0,200) # Residual standard deviation
    sigma_res_dsig[s] ~ dt(0, 0.01, 1)T(0,10) # Residual standard deviation
    sigma_res_asig[s] ~ dt(0, 0.01, 1)T(0,1) # Residual standard deviation
  
    tau_res_csig[s] <- 1/(sigma_res_csig[s]*sigma_res_csig[s])
    tau_res_dsig[s] <- 1/(sigma_res_dsig[s]*sigma_res_dsig[s])
    tau_res_asig[s] <- 1/(sigma_res_asig[s]*sigma_res_asig[s])
  
  }
  
  # csig site
    sigma_csig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_csig_site <- pow(sigma_csig_site, -2)

  # dsig site
    sigma_dsig_site ~ dt(0, 0.01, 1)T(0,20)
    tau_dsig_site <- pow(sigma_dsig_site, -2)

  # asig site
    sigma_asig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_asig_site <- pow(sigma_asig_site, -2)
  
  

}", fill=TRUE)
sink()


# Initial values ----------------------------------------------------------


# init_f <-  function(){
#   list(a = rnorm(data$N_an_sp_n,150,15),
#        d = rnorm(data$N_an_sp_n,0,1),
#        c = rnorm(data$N_an_sp_n,0,1),
#        
#        random_csig_site = rep(0,data$N_site),
#        sigma_csig_site = runif(1,0,100),
#        
#        random_asig_site = rep(0,data$N_site),
#        sigma_asig_site = runif(1,0,1),
#        
#        random_dsig_site = rep(0,data$N_site),
#        sigma_dsig_site = runif(1,0,15)
#        
#   )
# }

# inits <- list(init1 = init_f(),init2=init_f(), init3=init_f())



# Run models --------------------------------------------------------------

# 1
parameters <- c("asig","csig","dsig",
                
                "d","a","c",
                
                "sigma_csig_site",
                "sigma_dsig_site",
                "sigma_asig_site",
                "random_csig_site",
                "random_dsig_site",
                "random_asig_site", "deviance"
)

md_1 <- jags(data = data,
             parameters.to.save = parameters,
             model.file = "model_simul",
             # inits = inits,
             n.chains = 3,
             n.iter = 10000,
             n.burnin = 3000)

#saveRDS(md_1,"md_1_pheno.rds)


# library plot ------------------------------------------------------------

library(tidyverse)
library(bayesplot)

mat_md <- as.matrix(as.mcmc(md_1))


# csig --------------------------------------------------------------------


# jointure tableau pour identification espèce année
df_density <- mat_md%>%
  as.data.frame()%>%
  # select(starts_with("csig"))%>%
  select(starts_with("c["))%>%
  pivot_longer(everything(),
               names_to = "parameter", values_to = "est")%>%
  mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
  left_join(prod_f[,c("sitename", "year","species", "an_sp_n")])%>%
  distinct(est,an_sp_n,sitename,year,species, .keep_all = TRUE)

# plot des densités, un peu illisible
df_density%>%
  ggplot(aes(x = est,fill = parameter, group = parameter))+
  geom_density(alpha = 0.5)+
  facet_grid(species~year, 
             switch = "y")+
  theme_minimal() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )


# get mean value for csig (per species and per year)
df_density%>%
  group_by(year, species)%>%
  filter(year < 2020)%>%
  summarize(mean_csig = mean(est),
            sd__csig = sd(est))

# evolutin csig par année (plot)
df_density%>%
  group_by(year, species)%>%
  filter(year < 2020)%>%
  summarize(mean_csig = mean(est),
            sd__csig = sd(est))%>%
  ggplot(aes(x = year, y = mean_csig, color = species))+
  geom_point()+
  theme_bw()

# boxplot estimates per year and per species
df_density%>%
  ggplot(aes(x = year, y = est, fill = species, color = year))+
  geom_boxplot()+
  theme_bw()



# asig --------------------------------------------------------------------


# jointure tableau pour identification espèce année
df_density <- mat_md%>%
  as.data.frame()%>%
  # select(starts_with("csig"))%>%
  select(starts_with("asig["))%>%
  pivot_longer(everything(),
               names_to = "parameter", values_to = "est")%>%
  mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
  left_join(prod_f[,c("sitename", "year","species", "an_sp_n")])%>%
  distinct(est,an_sp_n,sitename,year,species, .keep_all = TRUE)

# plot des densités, un peu illisible
df_density%>%
  ggplot(aes(x = est,fill = parameter, group = parameter))+
  geom_density(alpha = 0.5)+
  facet_grid(species~year, 
             switch = "y")+
  theme_minimal() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )


# get mean value for csig (per species and per year)
df_density%>%
  group_by(year, species)%>%
  filter(year < 2020)%>%
  summarize(mean_asig = mean(est),
            sd__asig = sd(est))

# evolutin csig par année (plot)
df_density%>%
  group_by(year, species)%>%
  filter(year < 2020)%>%
  summarize(mean_asig = mean(est),
            sd__asig = sd(est))%>%
  ggplot(aes(x = year, y = mean_asig,fill = species, color = species))+
  geom_point()+
  theme_bw()+
  stat_smooth(method = 'lm')







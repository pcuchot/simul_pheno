# -----------------------------------------------------------------------------
# title : 1_simulate_data
# author : Paul Cuchot  
# date : 25/04/2023
# note : 
# -----------------------------------------------------------------------------


# Packages ----------------------------------------------------------------
library(tidyverse)


# -------------------------------------------------------------------------
# 1 Estimate phenology without taking into account young dispersal --------
# -------------------------------------------------------------------------

# arguments ? 

# n_session
# n_breeders
# n_years
# mean_ld (per year)

# 


# number of capture session per site
n_session <- 4

# number of breeding birds per site (catchable adults)
n_breeders <- 5

# number of breeding attempts
n_bredd_att <-  1

# breeding date (distrib to sample for breeding dates)
mean_ld <- 120
sd_ld <- 10

# sample n_breeders laying events
ld_dates <- round(rnorm(n_breeders, mean = mean_ld, sd = sd_ld))

# number of eggs per pair
mean_eggs <- 8
sd_eggs <- 5

n_eggs <- abs(round(rnorm(n_breeders, mean_eggs, sd_eggs)))

# fledglings dates 
fledgl_dates <- round(ld_dates)+35

# create a dataframe (one row per breeding pair)
df_breed <- data.frame(
  ld_dates = ld_dates,
  n_eggs = n_eggs,
  fledgl_dates = fledgl_dates
)

# Sample design -----------------------------------------------------------

# French CES = 4 capture sessions from day 120 to 180
t_capt <- round(seq(120,180, length.out = 4))

# Number of sample individual per session (~capture effort)
mean_n_capt <- 9

# Dataframe with n_adults and n_juveniles captured per session
df_session <- data.frame(t = t_capt,
                         n_capt_adults = NA,
                         n_capt_juveniles = NA) 

for(i in t_capt){
  # catchable adults
  n_adults <- n_breeders * 2
  # catchable juveniles
  n_juveniles <- sum(df_breed[df_breed$fledgl_dates < i ,]$n_eggs)
  
  capt_indiv <- sample( # how many individuals captured
    c(rep(0,n_adults), 
      rep(1,n_juveniles)),
    round(rnorm(1,mean_n_capt)), replace = TRUE)
  
  df_session[df_session$t ==i,"n_capt_adults"] <- sum(capt_indiv == 0)
  df_session[df_session$t ==i,"n_capt_juveniles"] <- sum(capt_indiv == 1) 
  
}

# calculate productivity per session
df_session <- df_session%>%
  mutate(prod = n_capt_juveniles/(n_capt_adults+n_capt_juveniles))

df_session


# Plot productivity -------------------------------------------------------

df_session%>%
  ggplot(aes(x = t, y = prod))+
  geom_point()+
  theme_light()


# function to simulated data for a signe site, N times --------------------

# 
simul_data <- function(n_breeders,
                       n_years, 
                       n_session){
  
  # mean laying date (change between years)
  mean_ld <- round(rnorm(n_years, 120, 15)) # (real pheno)
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
    
    # time for capture session
    t_capt <- round(seq(120,180, length.out = n_session))
    
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
      
      capt_indiv <- sample( # how many individuals captured
        c(rep(0,n_adults), 
          rep(1,n_juveniles)),
        round(rnorm(1,mean_n_capt)), replace = TRUE)
      
      df_session[df_session$t ==i,"n_capt_adults"] <- sum(capt_indiv == 0)
      df_session[df_session$t ==i,"n_capt_juveniles"] <- sum(capt_indiv == 1) 
      
    }
    
    df_session <- df_session%>%
      mutate(prod = n_capt_juveniles/(n_capt_adults+n_capt_juveniles))
    
    df_site <- rbind(df_site, df_session)
  }
  
  return(df_site)
  
}


data1 <- simul_data(n_breeders = 5,
                    n_years = 5, 
                    n_session = 4)

data1%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  theme_light()

# -----------------------------------------------------------------------------
# title : 1_simulate_data
# author : Paul Cuchot  
# date : 25/04/
# note : 
# -----------------------------------------------------------------------------


# Packages ----------------------------------------------------------------

library(tidyverse)


# -------------------------------------------------------------------------
# 1 Estimate phenology without taking into account young dispersal --------
# -------------------------------------------------------------------------


# number of capture session per site
n_session <- 4

# number of breeding birds per site (catchable adults)
n_breeders <- 5

# number of breeding attempts
n_bredd_att <-  1

# breeding date per site 
mean_ld <- 120
sd_ld <- 10

# sample n_breeders laying events
ld_dates <- round(rnorm(n_breeders, mean = mean_ld, sd = sd_ld))

# number of eggs per pair
mean_eggs <- 8
n_eggs <- rpois(n_breeders, mean_eggs)

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
                       rpois(1,mean_n_capt), replace = TRUE)
  
  df_session[df_session$t ==i,"n_capt_adults"] <- sum(capt_indiv == 0)
  df_session[df_session$t ==i,"n_capt_juveniles"] <- sum(capt_indiv == 1) 
 
}

# calculate productivity per session
df_session <- df_session%>%
  mutate(prod = n_capt_juveniles/(n_capt_adults+n_capt_juveniles))

df_session

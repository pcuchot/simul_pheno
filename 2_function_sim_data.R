
# -----------------------------------------------------------------------------
# title : 2_function_sim_data
# author : Paul Cuchot  
# date : 27/04/2023
# note : Create a function to simulate data for a single site 
# -----------------------------------------------------------------------------

# assumptions :

# - same number of breeders between years
# - mean number of eggs = 8
# - All chicks fledges and survive
# - Adults survive 
# - between years,  sessions are on the same days 
# - Birds do not migrate or immigrate
# - same capture probability juveniles / adults 


# Packages ----------------------------------------------------------------
require(tidyverse)


# function to simulated data for a signe site, N times --------------------


simul_data <- function(n_breeders,
                       n_years, 
                       n_session, 
                       start_ces = 80,
                       end_ces = 120,
                       mean_ld_site = 90){
  
  # mean laying date (among breeding individuals - change between years)
  mean_ld <- round(rnorm(n_years, mean_ld_site, 15)) # (real pheno)
  
  # variance in laying date for each year,
  sd_ld <- 7
  
  # can vary between years (then transform as a vector)
  # mean_var_ld <- 10
  # sd_ld <- rnorm(n_years, mean_var_ld, 15)
  
  
  # number of eggs per pair
  
  # (should depend on the laying date)
  
  mean_eggs <- 8
  sd_eggs <- 2
  
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
                            sd = sd_ld))# or sd_ld[k]
    
    # fledglings dates 
    fledgl_dates <- round(ld_dates)+50
    
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
      capt_indiv <- sample( 
        c(rep(0,n_adults), # adults
          rep(1,n_juveniles)), # juveniles
        round(rnorm(1,mean_n_capt, sd = 3)), # number of capture 
        replace = TRUE) 
      
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


# 
# sim data with
data1 <- simul_data(
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
# 
# 
# # plot productivity~day
# 
# data1$capt_sess%>%
#   ggplot(aes(x = t, y = prod, color = as.character(year)))+
#   geom_point()+
#   # add leaying dates
#   geom_vline(data = data1$mean_ld_year,
#              aes(xintercept = mean_ld,
#                  color = as.character(year)), alpha = 0.8)+
#   geom_line(alpha = 0.3)+
#   theme_light()








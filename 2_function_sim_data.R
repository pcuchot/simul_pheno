# -----------------------------------------------------------------------------
# title : 2_function_sim_data
# author : Paul Cuchot  
# date : 27/04/2023
# note : Create a function to simulate data per site 
# -----------------------------------------------------------------------------

# assumptions

# -
# -
# -
# -
# -
# -




# Packages ----------------------------------------------------------------
library(tidyverse)



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
  
  return(list(capt_sess = df_site, 
              mean_ld_year = mean_ld))
  
  
}

# sim data
data1 <- simul_data(n_breeders = 5,
                    n_years = 10, 
                    n_session = 12)


# plot productivity~day

data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  geom_line(alpha = 0.3)+
  theme_light()







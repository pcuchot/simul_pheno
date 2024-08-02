
# -----------------------------------------------------------------------------
# title : 2_function_sim_data
# author : Paul Cuchot  
# date : 27/04/2023
# note : Create a function to simulate data for a single site 
# -----------------------------------------------------------------------------

# assumptions :

# - same number of breeders between years
# - mean number of eggs = 8
# - number of fledgings depends on lyaing dates
# - Adults survive 
# - between years,  sessions are on the same days 
# - Birds do not migrate or immigrate
# - same capture probability juveniles / adults 


# Packages ----------------------------------------------------------------
require(tidyverse)
require(truncnorm)

# function to simulated data for a sigle site, N times --------------------


simul_data <- function(n_breeders = 10, # number of pair
                       n_years = 10, 
                       n_session = 50, 
                       start_ces = 100,
                       end_ces = 200,
                       mean_ld_site = 90,
                       selection_stre_ld = -0.003,
                       sd_ld = 7){
  
  
  # mean laying date (among breeding individuals - change between years)
  # mean_ld <- round(rnorm(n_years, mean_ld_site, 15)) # (real pheno)
  
  # for selection/variance explorations
  mean_ld <- 100
  
  # variance in laying date for each year,
  # sd_ld <- 7
  
  # can vary between years (then transform as a vector)
  # mean_var_ld <- 10
  # sd_ld <- rnorm(n_years, mean_var_ld, 15)
  
  
  # number of eggs per pair
  # (should depend on the laying date)
  mean_eggs <- 8
  # sd_eggs <- 2 # no need when sampling from poisson
  
  
  # final data_set
  df_site <- data.frame(t = NA,
                        n_capt_adults = NA,
                        n_capt_juveniles = NA, 
                        prod = NA,
                        year = NA)[0,] 
  
  
  for(k in 1:n_years){
    
    ## PHENOLOGY ##
    
    # sample n_breeders laying events
    ld_dates <- round(rnorm(n_breeders, 
                            mean = mean_ld[k], 
                            sd = sd_ld))# or sd_ld[k]
    
    # fledglings dates (50 = incubation time + rising)
    # imply that they all fledge at the same time 
    
    fledgl_dates <- round(ld_dates)+40
    
    # number of eggs per pair
    
    ## FECUNDITY ##
    
    n_eggs <- rpois(n_breeders, 
                    lambda = mean_eggs*exp(selection_stre_ld*ld_dates))
    
    # hist(rpois(1000,lambda = mean_eggs*exp(-0.003*ld_dates)), breaks = 100)
    
    
    # create a dataframe (one row per breeding pair)
    df_breed <- data.frame(
      ld_date = ld_dates,
      n_egg = n_eggs,
      fledgl_dates = fledgl_dates
    )
    
    
    #### sample (as CES design) #### 
    
    # choose days for capture session
    t_capt <- round(seq(start_ces, end_ces, 
                        length.out = n_session))
    
    # Number of sample individual per session (~capture effort)
    # for now: does not vary along the season
    
    mean_n_capt <- 30
    
    # Dataframe with n_adults and n_juveniles captured per session
    df_session <- data.frame(t = t_capt,
                             n_capt_adults = NA,
                             n_capt_juveniles = NA,
                             prod = NA,
                             year = as.character(k)) 
    
    for(i in t_capt){
      
      # catchable adults (no variation of survival during the season)
      n_adults <- n_breeders * 2
      
      # catchable juveniles (no variation of survival during the season)
      n_juveniles <- sum(df_breed[df_breed$fledgl_date < i ,]$n_egg)
      
      # sample birds among available individuals
      capt_indiv <- sample( 
        c(rep(0,n_adults), # adults
          rep(1,n_juveniles)), # juveniles
        rpois(1,mean_n_capt),
        
        replace = TRUE # allow recapture
      ) 
      
      # --> adults and juveniles have the same capture probability
      
      
      # how many adults
      df_session[df_session$t == i,"n_capt_adults"] <- sum(capt_indiv == 0)
      
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



# SELECTION ---------------------------------------------------------------


# selection = 0

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 90,
                    selection_stre_ld = 0,
                    sd_ld = 7)

data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()

# selection = -0.003

data1 <- simul_data(n_breeders = 10, # number of pair
           n_years = 1, 
           n_session = 50, 
           start_ces = 100,
           end_ces = 200,
           mean_ld_site = 90,
           selection_stre_ld = -0.003,
           sd_ld = 7)

data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()


# selection = -0.01

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 100,
                    selection_stre_ld = -0.01,
                    sd_ld = 7)

data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()


# VARIANCE LD -------------------------------------------------------------

# HIGHer variance in LD +
# no selection

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 100,
                    selection_stre_ld = 0,
                    sd_ld = 12)

high_var_no_sel <- data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()+ theme(legend.position = "none")+
  ggtitle("selection = 0, sd_layingdate = 12")

# HIGH variance in LD +
# strong selection

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 100,
                    selection_stre_ld = -0.01,
                    sd_ld = 12)

high_var_str_sel <- data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()+ theme(legend.position = "none")+
  ggtitle("selection = -0.01, sd_layingdate = 12")

# LOWer variance in LD +
# no selection

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 100,
                    selection_stre_ld = 0,
                    sd_ld = 5)

low_var_no_sel <- data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()+ theme(legend.position = "none")+
  ggtitle("selection = 0, sd_layingdate = 5")

# LOWer variance in LD +
# strong selection

data1 <- simul_data(n_breeders = 10, # number of pair
                    n_years = 1, 
                    n_session = 50, 
                    start_ces = 100,
                    end_ces = 200,
                    mean_ld_site = 100,
                    selection_stre_ld = -0.01,
                    sd_ld = 5)

low_var_str_sel <- data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  # add leaying dates
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld,
                 color = as.character(year)), alpha = 0.8)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()+ theme(legend.position = "none")+
  ggtitle("selection = -0.01, sd_layingdate = 5")

gridExtra::grid.arrange(high_var_str_sel,
                        low_var_str_sel,
                        high_var_no_sel,
                        low_var_no_sel)

library(brms)
library(tidyverse)

# function to simulate data 
  # - selection with optimum 
simul_data <- function(n_breeders = 1000, # number of pair
                       n_years = 10, 
                       n_session = 100, 
                       start_ces = 100,
                       end_ces = 200,
                       mean_ld_site = 90,
                       selection_stre_ld = -0.003,
                       sd_ld = 7,
                       fact_omega = 2,
                       # mean number of eggs per pair
                       mean_eggs = 8){
  
  
  # mean laying date (among breeding individuals - change between years)
  mean_ld <- round(rnorm(n_years, mean_ld_site, 15)) # (real pheno)
  
  # for selection/variance explorations
  # mean_ld <- 100
  
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
    
    # fledglings dates (40 = incubation time + rising)
    # imply that they all fledge at the same time 
    
    fledgl_dates <- ld_dates+40
    
    # number of eggs per pair
    
    
    ## FECUNDITY ##
  
    # optimum
    omega2 <-  fact_omega*sd_ld  # peak width
    
    shiftopt <- 20  # shift of mean phenotype from optimum = delay in repro
    
    opt <- mean_ld[k] - shiftopt
    
    
    #fitness function
    fitness <- exp(-(ld_dates - opt)^2/(2*omega2^2))
    
    n_eggs <- rpois(n_breeders, 
                    lambda = mean_eggs*fitness)
    
    # hist(n_eggs)
    
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
    
    mean_n_capt <- 200
    
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

# simulate data
data1 <- simul_data(n_breeders = 1000, # number of pair
                    n_years = 1, 
                    n_session = 150, 
                    start_ces = 50,
                    end_ces = 200,
                    mean_ld_site = 90,
                    selection_stre_ld = -0.003,
                    sd_ld = 7,
                    fact_omega = 2,
                    # mean number of eggs per pair
                    mean_eggs = 8)

# plot simulated data
data1$capt_sess%>%
  ggplot(aes(x = t, y = prod, color = as.character(year)))+
  geom_point()+
  geom_point(data = data1$mean_ld_year,
             aes(x = mean_ld, y = rep(0,1),
                 fill = as.character(year)), 
             shape=23, color="black", size=3)+
  geom_line(alpha = 0.3)+
  ylim(c(0,1))+
  theme_light()+ theme(legend.position = "none")

# fit non linear model*%>%
fit_loss <- brm(
  bf(
    prod ~ pinf/(1+exp((tm-t)/b)), 
    pinf ~ 1, tm ~ 1, b ~ 1,
    nl = TRUE),
  data = data1$capt_sess, family = gaussian(link = "identity"),
  prior = c(
    prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
    prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 130, 40))
    prior(normal(3, 5), nlpar = "b")
  ),
  control = list(adapt_delta = 0.9))

fit_loss

# plot model predictions
plot(conditional_effects(fit_loss), points = TRUE)

# transform model into matrix of simulation
md_df <- as.data.frame(as.matrix(as.mcmc(fit_loss)))

# estimated mean breeding time after selection (from luis doc)
# t_m - 40 - b log(1-pinf)
md_df <- md_df%>%
  mutate(ld_after_selection = 
           b_tm_Intercept-40-(b_b_Intercept*log(1-b_pinf_Intercept)))

# variance (from luis doc)
# (pi²b²)/3
md_df <- md_df%>%
  mutate(est_var_ld = (pi^2*b_b_Intercept^2)/3)
# mean estimated variance = 37
# and simulated laying date variance = 49 (BEFORE SELECTION)


# plot productivity data, with tm and estimated mean laying date
md_df%>%
  ggplot()+
  # estimated laying date after selection
  geom_vline(aes(xintercept = mean(md_df$ld_after_selection)), col = "blue")+
  # simulated laying date 
  geom_vline(data = data1$mean_ld_year,
             aes(xintercept = mean_ld), color = "red")+
  # simulated data
  geom_point(data = data1$capt_sess, aes(x = t, y = prod), 
             color = "black")+
  geom_text(
    x = 55,
    y = 0.5,
    label = 'estimated LD\n after selection',
    # size = 4, # font size
    show.legend = FALSE,
    color = 'blue'
  )+
  geom_text(
    x = 90,
    y = 0.5,
    label = 'simulated LD',
    # size = 4, # font size
    show.legend = FALSE,
    color = 'red'
  )+
  ylim(c(0,0.8))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(50, 200, by = 10))+
  ggtitle("estimated ld = t_m - 40 - b log(1-pinf)",
          sub = "selection with optimum: OMEGA 2*sd_ld, sd_ld = 7")

  
  
 
  
  
  
  
  
  
  
  
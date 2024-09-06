# -----------------------------------------------------------------------------
# title : 12_power_analysis
# author : Paul Cuchot  
# date : 12/08/2024
# note : simulate 3 times, with varying number of site (10, 100, 300)
# include site effects, but not year !
# -----------------------------------------------------------------------------

# load packages -----------------------------------------------------------
library(brms)
library(tidyverse)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# NO SELECTION ------------------------------------------------------------
# -------------------------------------------------------------------------
simul_data <- function(n_breeders = 1000, # number of pair
                       n_session = 150, 
                       start_ces = 50,
                       end_ces = 200,
                       # mean_ld_site = 90,
                       # selection_stre_ld = -0.003,
                       sd_ld = sd_,
                       mean_ld = 90,
                       fact_omega = omeg,
                       # mean number of eggs per pair
                       mean_eggs = 8, 
                       shiftopt = 10){
  
  # final data_set
  df_site <- data.frame(t = NA,
                        n_capt_adults = NA,
                        n_capt_juveniles = NA, 
                        prod = NA,
                        year = NA)[0,] 
  
  ## PHENOLOGY ##
  
  # sample n_breeders laying events
  ld_dates <- round(rnorm(n_breeders, 
                          mean = mean_ld, 
                          sd = sd_ld))# or sd_ld[k]
  
  # fledglings dates (40 = incubation time + rising)
  # imply that they all fledge at the same time 
  
  fledgl_dates <- ld_dates+40
  
  ## FECUNDITY ##
  
  # optimum
  omega2 <-  fact_omega*sd_ld  # peak width
  opt <- mean_ld - shiftopt
  
  #fitness function
  fitness <- exp(-(ld_dates - opt)^2/(2*omega2^2))
  
  n_eggs <- rpois(n_breeders, 
                  lambda = mean_eggs*fitness)
  
  
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
  
  
  # df with mean ld per year ("real breeding date")
  df_mean_ld <- data.frame(mean_ld = mean_ld)
  
  
  return(list(capt_sess = df_site, 
              mean_ld_year = df_mean_ld,
              R_et = mean(n_eggs),
              # mu_star = mean(fitness*ld_dates),
              mu_star = mean_ld-(shiftopt/((fact_omega^1)+1)),
              # var_star = var(fitness*ld_dates)))
              var_star = ((fact_omega^2)/((fact_omega^2)+1))*(sd_ld^2)))
  
}

data_md <- data.frame()

for(n_sites in c(10,100,300)){ # total number of sites
  
  for(n_breed in 5){ # number of breeders
    
    for(n_sess in 3:10){ # number of sessions
      
      data_int <- data.frame()
      
      for(n_site in 1:n_sites){ # 
        
        mean_ld_site <- rnorm(1, mean = 90, sd = 5)
        
        # extracted from code 10
        data_s <- simul_data(n_breeders = n_breed, # number of pair
                             n_session = n_sess, 
                             start_ces = 50,
                             end_ces = 200,
                             sd_ld = 5,
                             mean_ld = mean_ld_site,
                             fact_omega = 100,
                             # mean number of eggs per pair
                             mean_eggs = 10, 
                             shiftopt = 10)
        
        data_s$capt_sess$site <- n_site
        data_s$capt_sess$n_sess <- n_sess
        data_s$capt_sess$n_breeds <- n_breed
        data_s$capt_sess$tot_site <- n_sites
        data_s$capt_sess$mean_ld <- mean_ld_site
        
        data_md_s <- data_s$capt_sess
        
        data_int <- rbind(data_int, data_md_s)
        
      }
      
      # bind data
      
      # fit model
      fit <- brm(
        bf(
          prod ~ pinf/(1+exp((tm-t)/b)), 
          pinf ~ 1, tm ~ 1 + (1|site), b ~ 1 ,
          nl = TRUE),
        data = data_int, family = gaussian(link = "identity"),
        prior = c(
          prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
          prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 125, 40))
          prior(normal(3, 5), nlpar = "b")),
        control = list(adapt_delta = 0.9),
        chains = 4)
      
      # model fit to dataframe
      md_df <- as.data.frame(as.matrix(as.mcmc(fit)))
      
      # estimated mean laying date
      estim_pheno <- md_df$b_tm_Intercept - 40 - (md_df$b_b_Intercept*log(1-md_df$b_pinf_Intercept))
      data_md_s$est_pheno <- mean(estim_pheno)
      
      # estimated variance
      estim_var <- ((pi^2)*(md_df$b_b_Intercept^2))/3
      data_md_s$est_var <- mean(estim_var)
      
      data_md <- rbind(data_md, data_md_s)
      
    } # number of sessions
    
  } # number of breeders
  
} # total number of site 

saveRDS(data_md, "data_power_analysis.rds")



# -------------------------------------------------------------------------
data_md <- readRDS("data_power_analysis.rds")


# plot phenology accuracy
data_md %>% 
  # calculte accuracy soustraction of estimated and real ld
  # mutate(accuracy = abs(mean_ld-est_pheno)) %>%
  mutate(accuracy = abs(90 - est_pheno)) %>%
  distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(.~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote(delta[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())

# plot variance accuracy
data_md %>% 
  # calculte accuracy soustraction of estimated and real ld
  mutate(accuracy = abs(25-est_var)) %>%
  # filter(est_var<80)%>%
  distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(.~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote({sigma^2}[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())


# -------------------------------------------------------------------------
# WITH SELECTION ------------------------------------------------------------
# -------------------------------------------------------------------------
#OMEGA =2 

data_md2 <- data.frame()

for(n_sites in c(10,100,300)){ # total number of sites
  
  for(n_breed in 3:10){ # number of breeders
    
    for(n_sess in 3:12){ # number of sessions
      
      data_int <- data.frame()
      
      for(n_site in 1:n_sites){ # 
        
        mean_ld_site <- rnorm(1, mean = 90, sd = 5)
        
        # extracted from code 10
        data_s <- simul_data(n_breeders = n_breed, # number of pair
                             n_session = n_sess, 
                             start_ces = 50,
                             end_ces = 200,
                             sd_ld = 5,
                             mean_ld = mean_ld_site,
                             fact_omega = 2,
                             # mean number of eggs per pair
                             mean_eggs = 10, 
                             shiftopt = 10)
        
        data_s$capt_sess$site <- n_site
        data_s$capt_sess$n_sess <- n_sess
        data_s$capt_sess$n_breeds <- n_breed
        data_s$capt_sess$tot_site <- n_sites
        data_s$capt_sess$mean_ld <- mean_ld_site
        data_s$capt_sess$mu_star <- data_s$mu_star
        
        data_md_s <- data_s$capt_sess
        
        data_int <- rbind(data_int, data_md_s)
        
      }
      
      # bind data
      
      # fit model
      fit <- brm(
        bf(
          prod ~ pinf/(1+exp((tm-t)/b)), 
          pinf ~ 1, tm ~ 1 + (1|site), b ~ 1 ,
          nl = TRUE),
        data = data_int, family = gaussian(link = "identity"),
        prior = c(
          prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
          prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 125, 40))
          prior(normal(3, 5), nlpar = "b")),
        control = list(adapt_delta = 0.9),
        chains = 4)
      
      # model fit to dataframe
      md_df <- as.data.frame(as.matrix(as.mcmc(fit)))
      
      # estimated mean laying date
      estim_pheno <- md_df$b_tm_Intercept - 40 - (md_df$b_b_Intercept*log(1-md_df$b_pinf_Intercept))
      data_md_s$est_pheno <- mean(estim_pheno)
      
      # estimated variance
      estim_var <- ((pi^2)*(md_df$b_b_Intercept^2))/3
      data_md_s$est_var <- mean(estim_var)
      
      data_md2 <- rbind(data_md2, data_md_s)
      
    } # number of sessions
    
  } # number of breeders
  
} # total number of site 

saveRDS(data_md2, "data_power_analysis_sel.rds")


# -------------------------------------------------------------------------
data_md2 <- readRDS("data_power_analysis_sel.rds")


# plot phenology accuracy
data_md2 %>% 
  # calculte accuracy soustraction of estimated and real ld
  # mutate(accuracy = abs(mean_ld-est_pheno)) %>%
  mutate(accuracy = abs(90 - est_pheno)) %>%
  distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(.~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote(delta[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())

# plot variance accuracy
data_md2 %>% 
  # calculte accuracy soustraction of estimated and real ld
  mutate(accuracy = abs(25-est_var)) %>%
  # filter(est_var<80)%>%
  distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(.~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote({sigma^2}[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())


# -------------------------------------------------------------------------
data_md <- readRDS("data_power_analysis.rds")%>%
  select(tot_site, mean_ld, est_pheno, n_sess, n_breeds) %>% 
  distinct(tot_site, mean_ld, est_pheno, n_sess, n_breeds, .keep_all = T) %>% 
  mutate(sel = 100)
  

data_md2 <- readRDS("data_power_analysis_sel.rds")%>%
  select(tot_site, mean_ld, est_pheno, n_sess, n_breeds) %>% 
  distinct(tot_site, mean_ld, est_pheno, n_sess, n_breeds, .keep_all = T) %>% 
  mutate(sel = 2)

data_md3 <- rbind(data_md, data_md2)
  

data_md3 %>% 
  # calculte accuracy soustraction of estimated and real ld
  # mutate(accuracy = abs(mean_ld-est_pheno)) %>%
  mutate(accuracy = abs(90 - est_pheno)) %>%
  # distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(sel~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site),
                                     rows = tilde(omega)==.(sel)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote(delta[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())

# site level attempt

data_md3 %>% 
  # calculte accuracy soustraction of estimated and real ld
  # mutate(accuracy = abs(mean_ld-est_pheno)) %>%
  mutate(accuracy = abs(90 - est_pheno)) %>%
  # distinct(n_sess, n_breeds, tot_site, .keep_all = T) %>% 
  # filter(tot_site== 10) %>% 
  ggplot(aes(x = n_breeds, y = n_sess, fill = accuracy))+
  geom_tile()+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "olivedrab3",
                       high = "tomato1")+
  coord_fixed() +
  facet_grid(sel~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site),
                                     rows = tilde(omega)==.(sel)))+
  labs(y = "Number of capture session",
       x = "Number of breeding pairs", 
       fill = bquote(delta[simulated-estimated]))+
  theme_bw()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())

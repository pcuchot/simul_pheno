# -----------------------------------------------------------------------------
# title : 12bis_power_analysis
# author : Paul Cuchot  
# date : 04/09/2024
# note : update with brouillon on my desk
# -----------------------------------------------------------------------------

# load packages -----------------------------------------------------------
library(brms)
library(tidyverse)
# -------------------------------------------------------------------------


data_md <- data.frame()

n_sites = 100

# for(sd_betw in c(3,7)){
for(sd_betw in seq(1,10, length.out = 30)){
    
  data_int <- data.frame()
  
  for(n_site in 1:n_sites){ # 
    
    mean_ld_site <- rnorm(1, mean = 90, sd = sd_betw)
    
    # extracted from code 10
    data_s <- simul_data(n_breeders = 8, # number of pair
                         n_session = 10, 
                         start_ces = 50,
                         end_ces = 200,
                         sd_ld = 5,
                         mean_ld = mean_ld_site,
                         fact_omega = 4,
                         # mean number of eggs per pair
                         mean_eggs = 10, 
                         shiftopt = 10)
    
    data_s$capt_sess$site <- n_site
    data_s$capt_sess$sd_betw <- sd_betw
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
      prior(normal(3, 5), nlpar = "b"),
      prior(cauchy(5,2), class = "sd", nlpar = 'tm', group = "site")),
    control = list(adapt_delta = 0.9),
    chains = 4)
  
  # model fit to dataframe
  md_df <- as.data.frame(as.matrix(as.mcmc(fit)))
  
  # estimated mean laying date
  estim_pheno <- md_df$b_tm_Intercept - 40 - (md_df$b_b_Intercept*log(1-md_df$b_pinf_Intercept))
  data_md_s$est_pheno <- mean(estim_pheno)
  
  # mean between site random effect 
  data_md_s$rd_site <- as.numeric(mean(md_df$sd_site__tm_Intercept))
  
  # estimated variance
  data_md_s$est_var <- mean(((pi^2)*(md_df$b_b_Intercept^2))/3)
  
  data_md_s <- data_md_s %>% 
    distinct(est_pheno, .keep_all = T)
  
  # bind data
  data_md <- rbind(data_md, data_md_s)
  
}


saveRDS(data_md, "data_power_analysis2.rds")


# -------------------------------------------------------------------------
data_md <- readRDS("data_power_analysis2.rds")


# plot phenology accuracy
data_md %>% 
  ggplot(aes(x = sd_betw  , y = rd_site))+
  geom_point()+
  geom_function(fun = function(w){w},
                col = "grey67", size = 0.6, linetype = "solid")+
  # coord_fixed() +
  facet_grid(.~tot_site,
             labeller = label_bquote(cols = "N sites"==.(tot_site)))+
  labs(y = "Estimated between_site sd",
       x = "Simulated between_site sd")+
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

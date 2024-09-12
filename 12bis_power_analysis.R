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

n_sites = 150

for(n_sess in c(10,5)){ # UK, FRP
  
  for(sd_ld_ in c(4,8)){
    
    for(sd_betw in seq(0,10, length.out = 30)){
      
      data_int <- data.frame()
      
      for(n_site in 1:n_sites){ # 
        
        mean_ld_site <- rnorm(1, mean = 90, sd = sd_betw)
        
        # extracted from code 10
        data_s <- simul_data(n_breeders = 8, # number of pair
                             n_session = n_sess, 
                             start_ces = 50,
                             end_ces = 200,
                             sd_ld = sd_ld_,
                             mean_ld = mean_ld_site,
                             fact_omega = 4,
                             # mean number of eggs per pair
                             mean_eggs = 10, 
                             shiftopt = 10)
        
        data_s$capt_sess$site <- n_site
        data_s$capt_sess$sd_betw <- sd_betw
        data_s$capt_sess$tot_site <- n_sites
        data_s$capt_sess$mean_ld <- mean_ld_site
        data_s$capt_sess$sl_ld <- sd_ld_
        data_s$capt_sess$n_sess <- n_sess
        
        data_md_s <- data_s$capt_sess
        
        data_int <- rbind(data_int, data_md_s)
        
      }
      
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
          prior(cauchy(0.2,2), class = "sd", nlpar = 'tm', group = "site")),
        # control = list(adapt_delta = 0.9),
        control = list(adapt_delta = 0.95),
        chains = 4)
      
      # model fit to dataframe
      md_df <- as.data.frame(as.matrix(as.mcmc(fit)))
      
      # estimated mean laying date
      estim_pheno <- md_df$b_tm_Intercept - 40 - (md_df$b_b_Intercept*log(1-md_df$b_pinf_Intercept))
      data_md_s$est_pheno <- mean(estim_pheno)
      
      # mean between site random effect 
      data_md_s$mean_rd_site <- as.numeric(mean(md_df$sd_site__tm_Intercept))
      data_md_s$med_rd_site <- as.numeric(median(md_df$sd_site__tm_Intercept))
      
      # estimated variance
      data_md_s$est_var <- mean(((pi^2)*(md_df$b_b_Intercept^2))/3)
      data_md_s$sd_betw <- sd_betw
      
      data_md_s <- data_md_s %>% 
        distinct(est_pheno, .keep_all = T)%>%
        dplyr::select(tot_site, sl_ld, n_sess,sd_betw,mean_rd_site, med_rd_site, est_var, est_pheno)
      
      # bind data
      data_md <- rbind(data_md, data_md_s)
    }
  }
}

saveRDS(data_md, "data_power_analysis2.rds")


# -------------------------------------------------------------------------
data_md <- readRDS("data_power_analysis2.rds")


# plot variance 
var_plot <- data_md %>% 
  ggplot(aes(x = sd_betw, y = med_rd_site, color = as.factor(sl_ld)))+
  geom_line(size = 0.8)+
  geom_function(fun = function(w){w},
                col = "grey67", size = 0.6, linetype = "solid")+
  # coord_fixed() +
  facet_grid(.~n_sess,
             labeller = label_bquote(cols = "N sessions"==.(n_sess)),
             space = "fixed")+

  theme_bw()+
  labs(x = "Simulated between site sd",
       y = "Estimated between site sd",
       color = "Within site sd")+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  scale_color_viridis_d(end = 0.8)
  # coord_equal()
  # coord_fixed()


# plot grand mean
pheno_plot <- data_md %>% 
  ggplot(aes(x = sd_betw, y = est_pheno, color = as.factor(sl_ld)))+
  geom_line(size = 0.8)+
  geom_function(fun = function(w){90},
                col = "grey67", size = 0.6, linetype = "solid")+
  # coord_equal() +
  facet_grid(.~n_sess,
             labeller = label_bquote(cols = "N sessions"==.(n_sess)),
             space = "fixed")+
  
  theme_bw()+
  labs(x = "",
       y = "Estimated phenology",
       color = "Within site sd")+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "None")+
  scale_color_viridis_d(end = 0.8)

gridExtra::grid.arrange(pheno_plot, var_plot, 
                        nrow = 2)

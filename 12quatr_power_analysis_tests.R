# -----------------------------------------------------------------------------
# title : 12quatr_power_analysis_tests
# author : Paul Cuchot  
# date : 04/09/2024
# note : - random site for slope
#        - lot of capture sessions 
# -----------------------------------------------------------------------------

# load packages -----------------------------------------------------------
library(brms)
library(tidyverse)
# -------------------------------------------------------------------------


data_md <- data.frame()

n_sites = 100
rp <- 1 # number of repetition per simulation 

# for(n_sess in c(10,5)){ # UK, FRP
for(n_sess in c(70)){ # UK, FRP
  for(sd_ld_ in c(4,8)){ # within site variance in ld
  # for(sd_ld_ in c(4)){ # within site variance in ld
    
    for(sd_betw in seq(0,10, length.out = 10)){ # for a gradient of between site variance
      
      for (z in 1:rp){
        
        data_int <- data.frame()
        
        for(n_site in 1:n_sites){ # 
          
          
          mean_ld_site <- rnorm(1, mean = 90, sd = sd_betw)
          
          # extracted from code 10
          data_s <- simul_data(n_breeders = 8, # number of breeding pair per site
                               n_session = n_sess, 
                               start_ces = 50,
                               end_ces = 200,
                               sd_ld = sd_ld_,
                               mean_ld = mean_ld_site,
                               fact_omega = 100,
                               # mean number of eggs per pair
                               mean_eggs = 10, 
                               shiftopt = 10)
          
          data_s$capt_sess$site <- n_site
          data_s$capt_sess$sd_betw <- sd_betw
          data_s$capt_sess$tot_site <- n_sites
          data_s$capt_sess$mean_ld <- mean_ld_site
          data_s$capt_sess$sl_ld <- sd_ld_
          data_s$capt_sess$n_sess <- n_sess
          data_s$capt_sess$rep <- z
          
          data_md_s <- data_s$capt_sess
          
          data_int <- rbind(data_int, data_md_s)
          
        }
        
        # fit model
        fit <- brm(
          bf(
            prod ~ pinf/(1+exp((tm-t)/b)), 
            pinf ~ 1, tm ~ 1 + (1|site), b ~ 1,# + (1|site) ,
            nl = TRUE),
          data = data_int, family = gaussian(link = "identity"),
          prior = c(
            prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
            prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 125, 40))
            prior(normal(3, 5), nlpar = "b"),
            prior(cauchy(0.2,2), class = "sd", nlpar = 'tm', group = "site")),
          # control = list(adapt_delta = 0.9),
          control = list(adapt_delta = 0.9),
          chains = 3)
        
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
        
        # b
        data_md_s$mean_b <- as.numeric(mean(md_df$b_b_Intercept))
        data_md_s$med_b <- as.numeric(median(md_df$b_b_Intercept))
        
        # tm
        data_md_s$mean_tm <- as.numeric(mean(md_df$b_tm_Intercept))
        data_md_s$med_tm <- as.numeric(median(md_df$b_tm_Intercept))
        
        # pinf
        data_md_s$mean_pinf <- as.numeric(mean(md_df$b_pinf_Intercept))
        data_md_s$med_pinf <- as.numeric(median(md_df$b_pinf_Intercept))
        
        # make df
        data_md_s <- data_md_s %>% 
          distinct(est_pheno, .keep_all = T)%>%
          dplyr::select(tot_site, sl_ld, n_sess,sd_betw,mean_rd_site, 
                        med_rd_site, est_var, est_pheno,
                        mean_b, med_b, mean_tm, med_tm, 
                        mean_pinf, med_pinf) %>% 
          mutate(rp_ = z)
        
        # bind data
        data_md <- rbind(data_md, data_md_s)
        print(paste("n_sess =", n_sess))
        print(paste("sd_ld_ =", sd_ld_))
        print(paste("sd_betw =", sd_betw))
        print(paste("z =", z))
        
        # model with random slope
        # saveRDS(data_md, "data_power_analysis_rd_slope.rds")
        
        # model lot of sessions
        saveRDS(data_md, "data_power_analysis_lot_sess.rds")
        
      }
    }
  }
}



# -------------------------------------------------------------------------
data_md <- readRDS("data_power_analysis_rd_slope.rds")



# plot grand mean
pheno_plot <- data_md %>% 
  group_by(n_sess, sl_ld, sd_betw) %>%
  summarise(mean_est_pheno = mean(est_pheno),
            min_est_pheno = mean(est_pheno)- sd(est_pheno),
            max_est_pheno = mean(est_pheno)+ sd(est_pheno)) %>% 
  # mutate(rp2 = paste(n_sess, sl_ld, rp_)) %>% 
  ggplot(aes(x = sd_betw, #group = rp2,
             color = as.factor(sl_ld)))+
  
  geom_ribbon(aes(ymin =  min_est_pheno, 
                  ymax = max_est_pheno,
                  fill = as.factor(sl_ld)), alpha = 0.3) + 
  
  geom_line(aes(y = mean_est_pheno), size = 0.8)+
  geom_function(fun = function(w){90},
                col = "grey67", size = 0.6, linetype = "dashed")+
  # coord_equal() +
  facet_grid(.~n_sess,
             labeller = label_bquote(cols = "N sessions"==.(n_sess)),
             space = "fixed")+
  
  theme_bw()+
  labs(x = "",
       y = "Estimated phenology",
       color = "Within site sd",
       fill = "Within site sd")+
  theme(strip.placement = "outside",
        # strip.text = "",
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "None"
  )+
  scale_color_viridis_d(end = 0.8)+
  scale_fill_viridis_d(end = 0.8)

pheno_plot


# plot variance 

var_plot <- data_md %>% 
  group_by(n_sess, sl_ld, sd_betw) %>%
  summarise(mean_med_rd_site = mean(med_rd_site),
            min_med_rd_site = mean(med_rd_site)- sd(med_rd_site),
            max_med_rd_site = mean(med_rd_site)+sd(med_rd_site)) %>% 
  # mutate(rp2 = paste(n_sess, sl_ld, rp)) %>% 
  ggplot(aes(x = sd_betw, #group = rp2,
             color = as.factor(sl_ld)))+
  geom_ribbon(aes(ymin =  min_med_rd_site, 
                  ymax = max_med_rd_site,
                  fill = as.factor(sl_ld)), alpha = 0.3) + 
  geom_line(aes(, y = mean_med_rd_site), size = 0.8)+
  geom_function(fun = function(w){w},
                col = "grey67", size = 0.6, linetype = "dashed")+
  # coord_fixed() +
  facet_grid(.~n_sess,
             labeller = label_bquote(cols = "N sessions"==.(n_sess)),
             space = "fixed")+
  
  theme_bw()+
  labs(x = "",
       y = "Estimated between site sd",
       color = "Within site sd",
       fill = "Within site sd")+
  theme(strip.placement = "outside",
        # strip.text = "",
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "None"
  )+
  scale_color_viridis_d(end = 0.8)+
  scale_fill_viridis_d(end = 0.8)

var_plot



# intra site variance estimation 

intra_var_plot <- data_md %>% 
  group_by(n_sess, sl_ld, sd_betw) %>%
  mutate(est_var = sqrt(est_var)) %>% 
  summarise(mean_est_intv = mean(est_var),
            min_est_intv  = mean(est_var)- sd(est_var),
            max_est_intv  = mean(est_var)+ sd(est_var)) %>% 
  # mutate(rp2 = paste(n_sess, sl_ld, rp_)) %>% 
  ggplot(aes(x = sd_betw, #group = rp2,
             color = as.factor(sl_ld)))+
  
  geom_ribbon(aes(ymin =  min_est_intv, 
                  ymax = max_est_intv,
                  fill = as.factor(sl_ld)), alpha = 0.3) + 
  
  geom_line(aes(y = mean_est_intv), size = 0.8)+
  # geom_function(fun = function(w){4},
  #               col = "purple", size = 0.6, linetype = "dashed")+
  # geom_function(fun = function(w){8},
  #               col = "forestgreen", size = 0.6, linetype = "dashed")+
  # coord_equal() +
  facet_grid(.~n_sess,
             labeller = label_bquote(cols = "N sessions"==.(n_sess)),
             space = "fixed")+
  
  theme_bw()+
  labs(x = "Simulated between site sd",
       y = "Estimated within site sd",
       color = "Within site sd",
       fill = "Within site sd")+
  theme(strip.placement = "outside",
        # strip.text = "",
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "None"
  )+
  theme(strip.placement = "outside",
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  scale_color_viridis_d(end = 0.8)+
  scale_fill_viridis_d(end = 0.8)


intra_var_plot


gridExtra::grid.arrange(pheno_plot, var_plot, intra_var_plot,
                        nrow = 3)





# -----------------------------------------------------------------------------
# title : 10_explore_param_optim
# author : Paul Cuchot  
# date : 07/08/2024
# note : 
# -----------------------------------------------------------------------------


# load packages -----------------------------------------------------------
library(brms)
library(tidyverse)
# -------------------------------------------------------------------------


# function to simulate data 
# - selection with optimum 
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



# Assess phenology estimates with varying sigma and selection -------------
shiftopt = 20
k = 1
z = 1
df_rec <- list()
n_rep <- 10

for(omega_til in c(2,3,5, 100)){
  
  for(n_eggs in c(4,8,15)){ # mean laying date
    
    for(sd_ in seq(3,10, by = 0.1)){  # variance in laying date
      # for(sd_ in c(2,5,10)){  # variance in laying date
      
      for(q in 1:n_rep){
        
        # simulate data
        data_sim <- simul_data(n_breeders = 1000, # number of pairs
                               n_session = 150, 
                               start_ces = 50,
                               end_ces = 200,
                               mean_ld = 90,
                               # selection_stre_ld = -0.003,
                               sd_ld = sd_,
                               fact_omega = omega_til,
                               # mean number of eggs per pair
                               mean_eggs = n_eggs, 
                               shiftopt = shiftopt) 
        
        
        md_sim <- try(nls(prod ~ pinf/(1+exp((tm-t)/b)),
                          data = data_sim$capt_sess,
                          start = list(tm = 120, b = 3, pinf = 0.6)))
        
        # few model don't converge, here is a way to skip them 
        if(class(md_sim) == "try-error"){
          
          tm = NA
          b = NA
          pinf = NA
          
          # create data frame recording everything
          df_rec[[k]] <- data.frame(
            # get model parameters
            tm = NA,
            b = NA,
            pinf = NA,
            
            # tm - 40
            tm_m40 = NA,
            
            # mu bar * (mean breeding date from model parameters)
            mu_bar_star = NA,
            
            # sigma bar * (variance from model parameters)
            var_bar_star = NA,
            
            # mu* (from fitness function)
            mu_bar = NA,
            
            # sigma* (from fitness function)
            var_bar = NA,
            
            # from observed data (after selection ) 
            R_hat = NA, 
            
            # when taking selection
            R_hat_sel = NA,
            
            # add parameters
            sim_sd_sl = sd_,
            n_eggs = n_eggs,
            omega_tilde = omega_til, 
            sim = z,
            R_star = data_sim$R_et,
            var_star = data_sim$var_star,
            mu_star = data_sim$mu_star,
            mu = 90)
          
          k = k+1
          
        }else{
          
          tm = as.data.frame(coef(md_sim))[1,1]
          b = as.data.frame(coef(md_sim))[2,1]
          pinf = as.data.frame(coef(md_sim))[3,1]
          
          # create data frame recording everything
          df_rec[[k]] <- data.frame(
            # get model parameters
            tm = tm,
            b = b,
            pinf = pinf,
            
            # tm - 40
            tm_m40 = tm - 40,
            
            # mu bar *  equation 13
            mu_bar_star = tm - 40 - b*log(1-pinf),
            
            # sigma bar * equation 14
            var_bar_star = ((pi^2)*b^2)/3,
            
            # knowing selection parameters
            # mu bar  equation 19
            mu_bar =  (tm - 40 - b*log(1-pinf))+(shiftopt)/((omega_til^2)+1),
            
            # sigma bar  equation 19
            var_bar = (((omega_til^2)+1)/(omega_til^2))*(((pi^2)*b^2)/3),
            
            # from observed data (after selection ) 
            R_hat = (2*pinf)/(1-pinf), 
            
            # when taking selection
            # R_hat_sel = n_eggs*sqrt(((omega_til^2)+1)/(omega_til^2))*exp(((shiftopt^2)/((((omega_til^2)+1)/(omega_til^2))*(((pi^2)*b^2)/3)))/(2*(omega_til^2+1))),
            R_hat_sel = n_eggs*sqrt(((omega_til^2)+1)/(omega_til^2))*exp(((shiftopt^2/ ((((omega_til^2)+1)/(omega_til^2))*(((pi^2)*b^2)/3)))/((((omega_til^2)+1)/(omega_til^2))*(((pi^2)*b^2)/3)))/(2*(omega_til^2+1))),
            
            # add parameters
            sim_sd_sl = sd_,
            n_eggs = n_eggs,
            omega_tilde = omega_til, 
            sim = z,
            R_star = data_sim$R_et,
            var_star = data_sim$var_star,
            mu_star = data_sim$mu_star,
            mu = 90)
          
          k = k+1
          
        }
      }
      z = z+1
    }
  }
}


df_simul3 <- bind_rows(df_rec)
# saveRDS(df_simul3, "data_simul.rds")


# -------------------------------------------------------------------------
# try to include real post seletcion parameters
# -------------------------------------------------------------------------

df_simul3 <- readRDS("data_simul.rds")

# mu phenology

df_simul3 %>% 
  dplyr::select(mu,tm_m40, mu_bar_star, mu_bar,
                omega_tilde, sim, sim_sd_sl,n_eggs,mu_star) %>%
  pivot_longer( cols = c("tm_m40", "mu_bar_star", "mu_bar",
                         "mu", "mu_star")) %>%
  group_by(sim, name,sim_sd_sl, omega_tilde, n_eggs)%>%
  summarize(est_pheno = mean(value)) %>% as.data.frame() %>%
  mutate(est_real = case_when(name %in% c("mu","mu_star")~"real",
                              TRUE ~ "est"),
         # name = fct_relevel(name, c("mu","mu_star","tm_m40","mu_bar_star","mu_bar"))) %>%
         name = fct_relevel(name, c("mu","mu_star","mu_bar","tm_m40","mu_bar_star"))) %>%
  ggplot(aes(group = name))+
  geom_line(aes(x = sim_sd_sl^2, y = est_pheno, 
                color = name, linetype = name, size = name))+
  # ylim(80,91)+
  labs(x = bquote({sigma^2} [laydate~simulated]),
       y = bquote(mu[laydate~estimated]),
       color = bquote(mu), 
       linetype = bquote(mu),
       size = bquote(mu)) +
  theme_bw()+
  
  facet_grid(n_eggs~ omega_tilde, 
             labeller = label_bquote(cols = tilde(omega)==.(omega_tilde),
                                     rows = widehat(W[max])==.(n_eggs)))+
  scale_color_manual(
    labels = c(bquote(mu),
               bquote(mu~"*"),
               bquote(widehat(mu[z])),
               bquote(widehat(mu~"*"[z])),
               bquote(t[m]-T[f])),
    # breaks = c("mu", "mu_star", "tm_m40","mu_bar_star","mu_bar"),
    breaks = c("mu","mu_star","mu_bar","mu_bar_star","tm_m40"),
    values=c("brown", "orange2", 
             "darkslateblue", "olivedrab3", "tomato1"))+
  
  scale_linetype_manual(
    labels = c(bquote(mu),
               bquote(mu~"*"),
               bquote(widehat(mu[z])),
               bquote(widehat(mu~"*"[z])),
               bquote(t[m]-T[f])),
    breaks = c("mu", "mu_star", "tm_m40","mu_bar_star","mu_bar"),
    values = c("dashed","dashed","solid","solid","solid"))+
  
  scale_size_manual(
    labels = c(bquote(mu),
               bquote(mu~"*"),
               bquote(widehat(mu[z])),
               bquote(widehat(mu~"*"[z])),
               bquote(t[m]-T[f])),
    breaks = c("mu", "mu_star", "tm_m40","mu_bar_star","mu_bar"),
    values = c(0.6,0.6,0.9,0.9, 0.9))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())



# variance
df_simul3 %>%
  dplyr::select(var_star, var_bar, var_bar_star, omega_tilde, 
                sim, sim_sd_sl, n_eggs) %>%
  pivot_longer( cols = c("var_star","var_bar", "var_bar_star")) %>% 
  group_by(sim, name, sim_sd_sl, omega_tilde, n_eggs)%>% 
  summarize(est_var = mean(value)) %>% as.data.frame() %>% 
  mutate(name = fct_relevel(name, c("var_star","var_bar",
                                    "var_bar_star"))) %>%
  ggplot(aes(group = name))+
  geom_function(fun = function(w){w}, 
                col = "grey67", size = 0.6, linetype = "solid")+
  geom_line(aes(x = sim_sd_sl^2, y = est_var, 
                color = name, 
                linetype = name, 
                size = name), alpha = 1)+
  # ylim(0,120)+
  labs(x = bquote({sigma^2} [laydate~simulated]),
       y = bquote({sigma^2} [laydate~estimated]), 
       color = bquote({sigma^2}),
       linetype = bquote({sigma^2}),
       size = bquote({sigma^2}))+
  theme_bw()+
  scale_color_manual(
    labels = c(bquote(sigma^2~"*"),
               bquote(widehat({sigma^2})),
               bquote(widehat({sigma^2~"*"}))),
    breaks = c("var_star","var_bar", "var_bar_star"),
    values=c("orange2", "darkslateblue", "olivedrab3"))+
  
  scale_linetype_manual(
    labels = c(bquote(sigma^2~"*"),
               bquote(widehat({sigma^2})),
               bquote(widehat({sigma^2~"*"}))),
    breaks = c("var_star","var_bar", "var_bar_star"),
    values = c("dashed","solid","solid"))+
  
  scale_size_manual(
    labels = c(bquote(sigma^2~"*"),
               bquote(widehat({sigma^2})),
               bquote(widehat({sigma^2~"*"}))),
    breaks = c("var_star","var_bar", "var_bar_star"),
    values = c(0.6,0.9,0.9))+
  facet_grid(n_eggs~ omega_tilde, 
             labeller = label_bquote(cols = tilde(omega)==.(omega_tilde),
                                     rows = widehat(W[max])==.(n_eggs)))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())


# R
df_simul3 %>% 
  dplyr::select(R_hat, R_hat_sel, omega_tilde, sim, sim_sd_sl, n_eggs,
                R_star) %>%
  pivot_longer( cols = c("R_hat", "R_hat_sel","R_star")) %>%
  mutate(name = fct_relevel(name, c("R_star","R_hat",
                                    "R_hat_sel"))) %>%
  group_by(sim, name, sim_sd_sl, omega_tilde, n_eggs)%>% 
  summarize(est_R = mean(value)) %>%
  ggplot(aes(group = name))+
  # mean number of eggs before selection
  geom_hline(aes(yintercept = n_eggs),
             col = "black", size = 0.5, linetype = "dashed")+
  
  geom_line(aes(x = sim_sd_sl^2, y = est_R, 
                color = name,
                size = name,
                linetype = name), alpha = 1)+
  
  # facet_grid(.~as.factor(omega_tilde))+
  ylim(0,50)+
  # xlim(15 ,100)+
  labs(x = bquote({sigma^2} [laydate~simulated]),
       y = bquote(widehat(R) [estimated]),
       color = bquote(widehat(R)),
       size = bquote(widehat(R)),
       linetype = bquote(widehat(R)))+
  theme_bw()+
  scale_color_manual(
    labels = c(bquote(R^"*"),
               bquote(widehat(W[max])),
               bquote(widehat(R[p[infinity]]))),
    breaks = c("R_star","R_hat","R_hat_sel"),
    values=c("orange2", "darkslateblue", 
             "olivedrab3"))+
  
  scale_linetype_manual(
    labels = c(bquote(R^"*"),
               bquote(widehat(W[max])),
               bquote(widehat(R[p[infinity]]))),
    breaks  = c("R_star","R_hat","R_hat_sel"),
    values = c("dashed","solid","solid"))+
  
  scale_size_manual(
    labels = c(bquote(R^"*"),
               bquote(widehat(W[max])),
               bquote(widehat(R[p[infinity]]))),
    breaks  = c("R_star","R_hat","R_hat_sel"),
    values = c(0.6,0.9,0.9))+
  
  facet_grid(n_eggs~ omega_tilde, 
             labeller = label_bquote(cols = tilde(omega)==.(omega_tilde),
                                     rows = widehat(W[max])==.(n_eggs)))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.background = element_blank())
  





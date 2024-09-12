# -----------------------------------------------------------------------------
# title : 13_method_figure
# author : Paul Cuchot  
# date : 12/09/2024
# note : 
# -----------------------------------------------------------------------------

# load packages -----------------------------------------------------------
library(brms)
library(tidyverse)
# -------------------------------------------------------------------------


data_s <- simul_data(n_breeders = 10000, # number of pair
                     n_session = 1000, 
                     start_ces = 50,
                     end_ces = 170,
                     sd_ld = 5,
                     mean_ld = 90,
                     fact_omega = 10,
                     # mean number of eggs per pair
                     mean_eggs = 10, 
                     shiftopt = 10)

md_1 <- brm(
  bf(
    prod ~ pinf/(1+exp((tm-t)/b)), 
    pinf ~ 1, tm ~ 1, b ~ 1,
    nl = TRUE),
  data = data_s$capt_sess, family = gaussian(link = "identity"),
  prior = c(
    prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
    prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 130, 40))
    prior(normal(3, 5), nlpar = "b")
  ),
  control = list(adapt_delta = 0.9))

plot(conditional_effects(md_1), points = TRUE)
pinf = summary(md_1)$fixed$Estimate[1]
tm = summary(md_1)$fixed$Estimate[2]
b = summary(md_1)$fixed$Estimate[3]


data_s2 <- simul_data(n_breeders = 10000, # number of pair
                     n_session = 1000, 
                     start_ces = 50,
                     end_ces = 170,
                     sd_ld = 15,
                     mean_ld = 90,
                     fact_omega = 10,
                     # mean number of eggs per pair
                     mean_eggs = 10, 
                     shiftopt = 10)

md_2 <- brm(
  bf(
    prod ~ pinf/(1+exp((tm-t)/b)), 
    pinf ~ 1, tm ~ 1, b ~ 1,
    nl = TRUE),
  data = data_s2$capt_sess, family = gaussian(link = "identity"),
  prior = c(
    prior(normal(0.7, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
    prior(normal(125, 40), nlpar = "tm"), # hist(rnorm(1000, 130, 40))
    prior(normal(3, 5), nlpar = "b")
  ),
  control = list(adapt_delta = 0.9))

plot(conditional_effects(md_2), points = TRUE)
pinf2 = summary(md_2)$fixed$Estimate[1]
tm2 = summary(md_2)$fixed$Estimate[2]
b2 = summary(md_2)$fixed$Estimate[3]

ggplot()+
  # geom_point(data = data_s$capt_sess,
  #            aes(x = t, y = prod), color = "orange")+
  geom_function(fun = function(x) {predict(md_1, newdata = data.frame(t = x))},
    colour = "orange", size = 0.9)+

  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 5), 
                   y =..scaled..*0.9), 
               fill = "orange", 
               color = "orange", size = 0.8, linetype = "dashed",
               alpha = 0.1)+ 
  
  # geom_point(data = data_s2$capt_sess, 
  #            aes(x = t, y = prod), color = "tomato")+
    
    geom_function(
      fun = function(x) {predict(md_2, newdata = data.frame(t = x))},
      colour = "tomato", size = 0.9)+
  
  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 13), 
                   y =..scaled..*0.4), 
               fill = "tomato", 
               color = "tomato", size = 0.8, linetype = "dashed",
               alpha = 0.1)+
  theme_classic()+
  labs(y = "Productivity",
       x = "Date")+ 
  annotate("text", x=65, y=0.25, label= "sd = 15", color = "tomato")+
  annotate("text", x=75, y=0.75, label= "sd = 5", color = "orange")
  
  











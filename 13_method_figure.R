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

# different within site variance in laying date

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

dif_var <- ggplot()+
  geom_vline(xintercept = 130, linetype = "dashed", color = "grey56")+
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey56")+
  # geom_point(data = data_s$capt_sess,
  #            aes(x = t, y = prod), color = "orange")+
  geom_function(fun = function(x) {predict(md_1, newdata = data.frame(t = x))},
    colour = "orange", size = 0.9)+

  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 5), 
                   y =..scaled..*0.9), 
               fill = "orange", 
               color = "orange", size = 0.5, linetype = "solid",
               alpha = 0.1)+ 
  
  # geom_point(data = data_s2$capt_sess, 
  #            aes(x = t, y = prod), color = "tomato")+
    
    geom_function(
      fun = function(x) {predict(md_2, newdata = data.frame(t = x))},
      colour = "tomato", size = 0.9)+
  
  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 13), 
                   y =..scaled..*0.4), 
               fill = "tomato", 
               color = "tomato", size = 0.5, linetype = "solid",
               alpha = 0.1)+
  theme_classic()+
  xlim(55,155)+ #ylim(0,1)+
  labs(y = "",
       x = "")+ 
  annotate("text", x=65, y=0.25, label= "sd = 15", color = "tomato", size = 3.2)+
  annotate("text", x=75, y=0.75, label= "sd = 5", color = "orange", size = 3.2)+
  geom_segment(aes(x = 90, y = 0.92,
                   xend = 130, yend = 0.92), col = "grey56", size = 2,
               linetype = "solid",
               arrow = arrow(length = unit(0.3, "inches"),
                             type = "open", angle = 30))+
  annotate("text", x=110, y=0.85, label= "Fledging time", 
           color = "grey35", size = 3.2)+
  geom_point(aes(x = tm2, y = pinf2/2), color = "tomato", size = 2.5)+
  geom_point(aes(x = tm, y = pinf/2),color = "orange", size = 2.5)+
  annotate("text", x = 55, y = 0.95, label = "A)", color = 'black')
  # geom_segment(aes(x = 135, xend = tm2,
  #                  y = 0.6, yend = pinf2/2), size = 0.6,
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed", angle = 20))+
  # geom_segment(aes(x = 135, xend = tm,
  #                  y = 0.6, yend = pinf/2), size = 0.6,
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed", angle = 20))+
  # annotate("text", x=147, y=0.66, 
  #          label= "Estimated \n fledging peaks", color = "black", size = 3.2)
  # geom_segment(aes(x = 90, y = 0.38,
  #              xend = tm2-1.5, yend = 0.38), col = "tomato", size = 1,
  #              linetype = "dashed",
  #              arrow = arrow(length = unit(0.1, "inches"), ends = "both",
  #                            type = "closed", angle = 30))+
  # 
  # geom_segment(aes(x = 90, y = 0.42,
  #                  xend = tm, yend = 0.42), col = "orange", size = 1,
  #              linetype = "dashed", 
  #              arrow = arrow(length = unit(0.1, "inches"), ends = "both",
  #                            type = "closed", angle = 30))
  
dif_var  


# different in reproductive rate


data_s3 <- simul_data(n_breeders = 10000, # number of pair
                     n_session = 200, 
                     start_ces = 50,
                     end_ces = 170,
                     sd_ld = 5,
                     mean_ld = 90,
                     fact_omega = 10,
                     # mean number of eggs per pair
                     mean_eggs = 1.5, 
                     shiftopt = 10)

md_3 <- brm(
  bf(
    prod ~ pinf/(1+exp((tm-t)/b)), 
    pinf ~ 1, tm ~ 1, b ~ 1,
    nl = TRUE),
  data = data_s3$capt_sess, family = gaussian(link = "identity"),
  prior = c(
    prior(normal(0.45, 0.1), nlpar = "pinf"), # hist(rnorm(1000, 0.7, 0.1))
    prior(normal(125, 30), nlpar = "tm"), # hist(rnorm(1000, 130, 40))
    prior(normal(3, 5), nlpar = "b")
  ),
  control = list(adapt_delta = 0.9))

plot(conditional_effects(md_3), points = TRUE)

pinf3 = summary(md_3)$fixed$Estimate[1]
tm3 = summary(md_3)$fixed$Estimate[2]
b3 = summary(md_3)$fixed$Estimate[3]


dif_prod <- ggplot()+
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey56")+
  geom_vline(xintercept = 130, linetype = "dashed", color = "grey56")+
  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 5), 
                   y =..scaled..*0.9), 
               fill = "orange", 
               color = "orange", size = 0.5, linetype = "solid",
               alpha = 0.1)+ 

  geom_density(aes(x = rnorm(1000000, mean = 90, sd = 5), 
                   y =..scaled..*0.4), 
               fill = "tomato", 
               color = "tomato", size = 0.5, linetype = "solid",
               alpha = 0.1)+
  
  geom_function(fun = function(x) {predict(md_3, newdata = data.frame(t = x))},
    colour = "tomato", size = 0.9)+
  
  geom_function(fun = function(x) {predict(md_1, newdata = data.frame(t = x))},
                colour = "orange", size = 0.9)+
  
  theme_classic()+
  theme(axis.text.y=element_blank())+
  xlim(55,155)+ #ylim(0,1)+
  labs(y = "",
       x = "")+ 
  annotate("text", x=65, y=0.25, label= "Max eggs = 1.5", 
           color = "tomato", size = 3.2)+
  annotate("text", x=70, y=0.75, label= "Max eggs = 10", 
           color = "orange", size = 3.2)+
  geom_segment(aes(x = 90, y = 0.92,
                   xend = 130, yend = 0.92), col = "grey56", size = 2,
               linetype = "solid",
               arrow = arrow(length = unit(0.3, "inches"),
                             type = "open", angle = 30))+
  annotate("text", x=110, y=0.85, label= "Fledging time", 
           color = "grey35", size = 3.2)+
  geom_point(aes(x = tm3, y = pinf3/2), 
             color = "tomato", size = 2.5)+
  geom_point(aes(x = tm, y = pinf/2), 
             color = "orange", size = 2.5)+
  annotate("text", x = 55, y = 0.95, label = "B)", color = 'black')
  # geom_segment(aes(x = 135, xend = tm3,
  #                  y = 0.6, yend = pinf3/2), size = 0.6,
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed", angle = 20))+
  # geom_segment(aes(x = 135, xend = tm,
  #                  y = 0.6, yend = pinf/2), size = 0.6,
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed", angle = 20))+
  # annotate("text", x=147, y=0.66, 
  #          label= "Estimated \n fledging peaks", 
  #          color = "black", size = 3.2)
  # geom_segment(aes(x = 90, y = 0.22,
  #                  xend = tm3, yend = 0.22), col = "tomato", size = 1,
  #              linetype = "dashed",
  #              arrow = arrow(length = unit(0.1, "inches"), ends = "both",
  #                            type = "closed", angle = 30))+
  # 
  # geom_segment(aes(x = 90, y = 0.40,
  #                  xend = tm, yend = 0.40), col = "orange", size = 1,
  #              linetype = "dashed", 
  #              arrow = arrow(length = unit(0.1, "inches"), ends = "both",
  #                            type = "closed", angle = 30))
dif_prod

gridExtra::grid.arrange(dif_var,dif_prod, nrow = 1,
                        bottom = 'Time in days', left = "Productivity")

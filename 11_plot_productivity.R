# -----------------------------------------------------------------------------
# title : 11_plot_productivity
# author : Paul Cuchot  
# date : 08/08/2024
# note : 
# -----------------------------------------------------------------------------


# plot simulated data -----------------------------------------------------


# productivity
ggplot()+
  # laying date after selection DISTRIBUTION
  geom_density(aes(x = rnorm(100000,
                             tm - 40 - b*log(1-pinf), 
                             ((pi^2)*b^2)/3),  
                   y = ..density..*4), col = "white", fill = "grey69", alpha = 0.5)+
  geom_point(data = data_sim$capt_sess,
             aes(x = t, y = prod, color = year_f))+
  geom_line(data = data_sim$capt_sess,
            aes(x = t, y = prod, color = year_f), alpha = 0.3)+
  # simulated mean laying date
  geom_point(data = data_sim$mean_ld_year,
             aes(x = mean_ld, y = rep(0,1)), fill = "blue", 
             shape=23, color="black", size=3)+
  # plot model
  geom_function(fun = sigmoid_f, color = "red")+
  
  # post selection distribution (back calculated from models)
  # tm - 40 - b*log(1-pinf)
  geom_point(data = data.frame(x = tm - 40 - b*log(1-pinf),
                               y = 0),
             aes(x = x, y = y), 
             shape=23, color="black", fill = 'red', size=3)+
  
  
  #laying date before selection (from known parameters)
  # geom_density(aes(x = rnorm(100000, data_sim$mean_ld_year$mean_ld, sd_),
  #                  y=..density..*3))+
  
  # laying date after selection (equation 16) from known selection parameters
  geom_density(aes(x = rnorm(100000,
                             data_sim$mean_ld_year$mean_ld - ((10)/((omeg)^2+1)),
                             (sd_^2)*((omeg*sd_)^2)/(((omeg*sd_)^2)+(sd_)^2)),
                   y = ..density..*4))+
  
  # 
  # ylim(c(0,0.9))+
  theme_light()+
  theme(legend.position = "none")

###





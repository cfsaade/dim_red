library(ggplot2)
library(reshape2)
data = read.csv('output/0_All_data.csv', header = T)


# computing the relative absolute error of mean, sd and acf_1
data$mean_RAE = abs(data$x_reduced_mean - data$R_mean)/pmax(data$R_mean, data$x_reduced_mean)
data$sd_RAE = abs(data$x_reduced_sd - data$R_sd)/pmax(data$R_sd, data$x_reduced_sd)
data$acf_1_RAE = abs(data$x_reduced_acf_1 - data$R_acf_1)/pmax(abs(data$R_acf_1), abs(data$x_reduced_acf_1))



ggplot(data = data, mapping = aes(x = noise_amplitude, y = x_reduced_sd, color = as.factor(noise_cov_dist))) +
  geom_point(alpha = 0.2) +
  facet_grid(~kind) +
  stat_summary(mapping = aes(x = noise_amplitude, y = R_sd, fill = as.factor(noise_cov_dist)),
               color = 'transparent',
               fun.min = function(z) { quantile(z,0) },
               fun.max = function(z) { quantile(z,1) },
               fun = median,
               geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = R_sd, color = as.factor(noise_cov_dist)),
    fun = median,
    geom = 'line')  +
  xlab('Noise amplitude') +
  ylab('Biomass standard deviation') +
  scale_color_discrete("Spatial autocorrelation of noise") +
  scale_fill_discrete("Spatial autocorrelation of noise") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave('sd_prediction.pdf', width = 210, height = 90, unit = 'mm')





ggplot(data = data, mapping = aes(x = noise_amplitude, y = x_reduced_sd, color = kind)) +
  geom_point(alpha = 0.2) +
  facet_grid(~noise_cov_dist) +
  stat_summary(mapping = aes(x = noise_amplitude, y = R_sd, fill = kind),
               color = 'transparent',
               fun.min = function(z) { quantile(z,0) },
               fun.max = function(z) { quantile(z,1) },
               fun = median,
               geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = R_sd, color = kind),
    fun = median,
    geom = 'line')  +
  xlab('Noise amplitude') +
  ylab('Biomass standard deviation') +
  scale_color_discrete("Spatial autocorrelation of noise") +
  scale_fill_discrete("Spatial autocorrelation of noise") +
  theme_bw() +
  theme(legend.position = "bottom")

## not used : ################################

ggplot(data = data, mapping = aes(x = noise_amplitude, y = R_sd, color = as.factor(noise_cov_dist))) +
  geom_point(alpha = 0.2) +
  facet_grid(~kind) +
  stat_summary(mapping = aes(x = noise_amplitude, y = R_sd, fill = as.factor(noise_cov_dist)),
               color = 'transparent',
               fun.min = function(z) { quantile(z,0) },
               fun.max = function(z) { quantile(z,1) },
               fun = median,
               geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = R_sd, color = as.factor(noise_cov_dist)),
    fun = median,
    geom = 'line') +
  xlab('Noise amplitude') +
  ylab('Biomass standard deviation') +
  scale_color_discrete("Spatial autocorrelation of noise") +
  scale_fill_discrete("Spatial autocorrelation of noise") +
  theme_bw() +
  theme(legend.position = "bottom")




ggplot(data = data, mapping = aes(x = noise_amplitude, y = x_reduced_sd, color = as.factor(noise_cov_dist))) +
  geom_point(alpha = 0.2) +
  facet_grid(~kind) +
  stat_summary(mapping = aes(x = noise_amplitude, y = x_reduced_sd, fill = as.factor(noise_cov_dist)),
               color = 'transparent',
               fun.min = function(z) { quantile(z,0) },
               fun.max = function(z) { quantile(z,1) },
               fun = median,
               geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = x_reduced_sd, color = as.factor(noise_cov_dist)),
    fun = median,
    geom = 'line') +
  xlab('Noise amplitude') +
  ylab('Biomass standard deviation') +
  scale_color_discrete("Spatial autocorrelation of noise") +
  scale_fill_discrete("Spatial autocorrelation of noise") +
  theme_bw() +
  theme(legend.position = "bottom")




ggplot(data = data, mapping = aes(x = noise_amplitude, y = x_reduced_sd, color = as.factor(noise_cov_dist))) +
  geom_point(alpha = 0.2) +
  facet_grid(~kind) +
  stat_summary(mapping = aes(x = noise_amplitude, y = R_sd, fill = as.factor(noise_cov_dist)),
               color = 'transparent',
               fun.min = function(z) { quantile(z,0) },
               fun.max = function(z) { quantile(z,1) },
               fun = median,
               geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = R_sd, color = as.factor(noise_cov_dist)),
    fun = median,
    geom = 'line')  +
  xlab('Noise amplitude') +
  ylab('Biomass standard deviation') +
  scale_color_discrete("Spatial autocorrelation of noise") +
  scale_fill_discrete("Spatial autocorrelation of noise") +
  theme_bw() +
  theme(legend.position = "bottom")


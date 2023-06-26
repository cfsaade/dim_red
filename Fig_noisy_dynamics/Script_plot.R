library(ggplot2)
library(reshape2)
data = read.csv('output/0_All_data.csv', header = T)


# computing the relative absolute error of mean, sd and acf_1
data$mean_RAE = abs(data$x_reduced_mean - data$R_mean)/pmax(data$R_mean, data$x_reduced_mean)
data$sd_RAE = abs(data$x_reduced_sd - data$R_sd)/pmax(data$R_sd, data$x_reduced_sd)
data$acf_1_RAE = abs(data$x_reduced_acf_1 - data$R_acf_1)/pmax(abs(data$R_acf_1), abs(data$x_reduced_acf_1))

data_melted = melt(data,
         measure.vars = c("mean_RAE", "sd_RAE", "acf_1_RAE"),
         variable.name = 'Metric',
         value.name = 'Relative_absolute_error')

metric.lab = c('Mean', 'Standard deviation', 'Autocorrelation lag 1')
names(metric.lab) = c('mean_RAE', 'sd_RAE', 'acf_1_RAE')


ggplot(data = data_melted, mapping = aes(x = noise_amplitude, y = Relative_absolute_error, color = as.factor(noise_cov_dist), shape = as.factor(noise_cov_dist))) +
  geom_point(alpha = 0.5) +
  facet_grid(Metric ~ kind, labeller = labeller(Metric = metric.lab)) + 
  stat_summary(
    mapping = aes(x = noise_amplitude, y = Relative_absolute_error, fill = as.factor(noise_cov_dist)),
    color = 'transparent',
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = noise_amplitude, y = Relative_absolute_error, color = as.factor(noise_cov_dist)),
    fun = median,
    geom = 'line') +
  scale_color_discrete('Spatial autocorrelation of noise') +
  scale_fill_discrete('Spatial autocorrelation of noise') +
  scale_shape_discrete('Spatial autocorrelation of noise') +
  xlab('Noise amplitude') +
  ylab('Relative absoulte error') +
  theme_bw() +
  theme(legend.position = 'bottom')

ggsave('RAE_noisy.pdf', width = 210, height = 148, unit = 'mm')




# some verification of the bias by plotting predicted vs. observed
ggplot(data = data, mapping = aes(x = R_mean, y = x_reduced_mean, shape = kind, color = kind)) +
  geom_point(alpha = 0.5) +
  xlim(c(25, 30)) +
  ylim(c(25, 30)) +
  theme_bw()

ggplot(data = data, mapping = aes(x = R_sd, y = x_reduced_sd, shape = kind, color = kind)) +
  geom_point(alpha = 0.5) +
  xlim(c(0, 6)) +
  ylim(c(0, 6)) +
  theme_bw()

ggplot(data = data, mapping = aes(x = R_acf_1, y = x_reduced_acf_1, shape = kind, color = kind)) +
  geom_point(alpha = 0.5) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  theme_bw()

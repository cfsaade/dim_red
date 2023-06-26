library(ggplot2)
data = read.csv('output/0_All_data.csv', header = T)



ggplot(data = data, mapping = aes(x = mu, y = abs(x_reduced - R)/max(R, x_reduced), shape = as.factor(n_patch), color = as.factor(n_patch))) +
  geom_point(alpha = 0.2) +
  stat_summary(
    mapping = aes(x = mu, y = abs(x_reduced - R)/max(R, x_reduced), fill = as.factor(n_patch)),
    color = 'transparent',
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom = 'ribbon', alpha = 0.3) +
  stat_summary(
    mapping = aes(x = mu, y = abs(x_reduced - R)/max(R, x_reduced), color = as.factor(n_patch)),
    fun = median,
    geom = 'line') +
  facet_grid(k_sd~kind) +
  scale_color_discrete('Network size') +
  scale_fill_discrete('Network size') +
  scale_shape_discrete('Network size') +
  xlab(expression(paste('Dispersal rate (', mu, ')'))) +
  ylab('Relative absoulte error') +
  theme_bw() +
  theme(legend.position = 'bottom')

ggsave('RAE_heterogeneity.pdf', width = 210, height = 148, unit = 'mm')

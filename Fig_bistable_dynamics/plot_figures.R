rm(list = ls())

library(ggplot2)
library(reshape2)
data = read.csv('output/0_All_data.csv', header = T)


# computing the relative absolute error of mean, sd and acf_1
#data$dist_pred = abs(data$x_reduced - data$R_predicted)/pmax(abs(data$R_predicted), abs(data$x_reduced)) # distance to predicted equilibrium
#data$dist_R_alt_1 = abs(data$x_reduced - data$R_alternative_1)/pmax(abs(data$R_alternative_1), abs(data$x_reduced)) # distance to eq 1
#data$dist_R_alt_2 = abs(data$x_reduced - data$R_alternative_2)/pmax(abs(data$R_alternative_2), abs(data$x_reduced)) # distance to eq 2
#data$dist_any = pmax(data$dist_R_alt_1, data$dist_R_alt_2)


## Defining a relative absolute functions between a pair of values that deals with values close to 0: if both values are < 10^-3, return 0
RAE = function(pair_values){
  if (max(abs(pair_values)) < 10^-3){
    return(0)
  }
  return(abs(pair_values[1] - pair_values[2])/max(abs(pair_values[1]), abs(pair_values[2])))
}


# computing RAEs
data$dist_pred = apply(data[,c("x_reduced", "R_predicted")], FUN = RAE, MARGIN = 1) # distance to predicted equilibrium
data$dist_R_alt_1 = apply(data[,c("x_reduced", "R_alternative_1")], FUN = RAE, MARGIN = 1) # distance to eq 1
data$dist_R_alt_2 = apply(data[,c("x_reduced", "R_alternative_2")], FUN = RAE, MARGIN = 1) # distance to eq 2
data$dist_any = pmin(data$dist_R_alt_1, data$dist_R_alt_2)

data_melted = melt(data,
         measure.vars = c("dist_pred", 'dist_any'),
         variable.name = 'Metric',
         value.name = 'Relative_absolute_error')

# label for panels
metric.lab = c('Distance to predicted eq.', 'distance to closest eq.')
names(metric.lab) = c('dist_pred', 'dist_any')


ggplot(data = data_melted, mapping = aes(x = mu, y = Relative_absolute_error, color = as.factor(c), shape = as.factor(c))) +
  geom_point(alpha = 0.3) +
  facet_grid(Metric ~ kind, labeller = labeller(Metric = metric.lab)) +
  stat_summary(
    mapping = aes(x = mu, y = Relative_absolute_error, fill = as.factor(c)),
    color = 'transparent',
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom = 'ribbon', alpha = 0.2) +
  stat_summary(
    mapping = aes(x = mu, y = Relative_absolute_error, color = as.factor(c)),
    fun = median,
    geom = 'line') +
  xlab(expression(paste('Dispersal rate (', mu, ')'))) +
  ylab('Relative absoulte error') +
  scale_color_discrete('Allee effect threshold') +
  scale_fill_discrete('Allee effect threshold') +
  scale_shape_discrete('Allee effect threshold') +
  theme_bw() +
  theme(legend.position = 'bottom')


ggsave('RAE_bistable.pdf', width = 210, height = 148, unit = 'mm')


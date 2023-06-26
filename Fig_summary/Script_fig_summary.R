# a script that tests the one-dimensional spectral reduction over a list of networks

# sources and libraries ####################################################################
rm(list = ls())
library(deSolve)
library(igraph)

# functions for the reduction
source("../one_dimension_spectral_reduction.R")

# functions for dynamics
source("../Dynamics/metapop_dynamics.R")
source("../Dynamics/Growth_functions.R")
source("../Dynamics/Dispersal_functions.R")
source("../Networks/Network_functions.R")

## declaring a function to plot a network
plot_network = function(network, size = 1){
  # extracting the number of patches
  n_patch = network$summary$n_patch
  
  # dataframe of point coordinates
  points = network$points_coords
  points$size = size**2
  
  # making a dataframe of link coordinates
  links = data.frame(x = 0, y = 0, link = 0, weight = 0)
  link = 0
  W = network$dispersal_matrix
  for(k in 1:(n_patch-1)){
    for(l in (k+1):n_patch){
      if(W[k,l]>0){
        link = link+1
        temp_data = data.frame(x = points[c(k, l), 1], y = points[c(k,l), 2], link = link, weight = (W[k,l] + W[l,k])/2)
        links = rbind(links, temp_data)
      }
    }
  }
  
  return(ggplot() +
           geom_line(data = links, mapping = aes(x = x, y = y, group = link), color = "grey", size = 0.5) +
           geom_point(data = points, mapping = aes(x = x, y = y, size = size), color = 'black') +
           labs(x="",y="",color="Biomass", size="Weight") +
           theme_classic() +
           labs(x="",y="") +
           theme(axis.text.x = element_blank(),legend.position="bottom",
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line = element_line(colour = "white")) +
           xlim(min(points$x) - 0.1, max(points$x) + 0.1) +
           ylim(min(points$y) - 0.1, max(points$y) + 0.1) +
           guides(size = "none") +
           xlim(0, 1) +
           ylim(0, 1) +
           scale_size(range = c(0, 3)))
}

## parameters for the model ###########################################
r = 1
K = 30
c = 0.3

par = list(r = r,
           K = K,
           c = c)


mu = 0.1

# timestep for the stochastic integration
dt = 0.01

## loading network ####################################################
network_id = 5
network = read_network(paste0("../Fig_noisy_dynamics/Network_list/RGG_", network_id))
disp_matrix = network$dispersal_matrix
W = disp_matrix*mu # scaling the disp_matrix by mu
n_patch = ncol(disp_matrix)
#n_patch = ncol(disp_matrix)
#disp_matrix = as.matrix(disp_matrix, nrow = n_patch)

g = graph_from_adjacency_matrix(network$adjacency_matrix)

# adjusting layout ####################################################
V(g)$label = ""
V(g)$size = 5
E(g)$arrow.size = 0
coords = layout.lgl(g)
plot(g, layout = coords)

plot_network(network)


# simulating the full dynamics #######################################

# generating inital density using a lognormal distribution of mean m = m and sd = s = 0.1
#m = 0.7
#s = 0.3
#location <- log(m^2 / sqrt(s^2 + m^2))
#shape <- sqrt(log(1 + (s^2 / m^2)))

#x_init = rlnorm(n_patch, location, shape)*K

x_init = runif(n_patch, 0, K)


# setting the time sequence
t = seq(1:1000)

# computing the odsr
odsr.object = odsr(W, type = 'stochastic', cov_matrix = network$cov_matrix) # computing the one dimensional spectral reduction

if (Re(odsr.object$alpha) < 0 ){
  odsr.object = odsr(W, type = 'stochastic', cov_matrix = network$cov_matrix, n_eig = 2)
}

odsr.object$a = sapply(odsr.object$a, Re)
odsr.object$alpha = Re(odsr.object$alpha)
odsr.object$alpha.out = Re(odsr.object$alpha.out)
odsr.object$beta = Re(odsr.object$beta)
odsr.object$cov = Re(odsr.object$cov)


# running the full simulation
out = ode(y = x_init,
          times = t,
          func = metapop_dynamics_stochastic,
          parms = par,
          growth_function = logistic_allee,
          dispersal_function = linear_dispersal,
          W = W,
          
          Cov = network$cov_matrix/dt,
          method = "euler",
          hini = 0.01)

eq_density = out[nrow(out), 2:ncol(out)]


# making a dataframe containing the points coordinates
data_points = data.frame(x = coords[,1], y = coords[,2], x_init = x_init, x_eq = eq_density)

# making a dataframe containing the links
data_links = data.frame(x = 0, y = 0, link = 0, weight = 0)
link = 0
for(k in 1:(n_patch-1)){
  for(l in (k+1):n_patch){
    if(W[k,l]>0){
      link = link+1
      temp_data = data.frame(x = coords[c(k, l), 1], y = coords[c(k,l), 2], link = link, weight = (W[k,l] + W[l,k])/2)
      data_links = rbind(data_links, temp_data)
    }
  }
}

data_links = data_links[-1,]


# panel a ##################################################################################################
# plotting initial state
panel_a = ggplot() +
  geom_line(data = data_links, mapping = aes(x = x, y = y, group = link, size = weight), color = "grey") +
  geom_point(data = data_points, mapping = aes(x = x, y = y, color = x_init), size = 5) +
  labs(x="",y="",color="Biomass", size="Weight") +
  scale_size_continuous(range = c(0.05,1)) +
  scale_colour_gradient(low = "red", high = "green", limits = c(0, K*1.2)) +
  theme_classic() +
  labs(x="",y="") +
  theme(axis.text.x = element_blank(),legend.position="bottom",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "white")) +
  xlim(min(data_points$x) - 1, max(data_points$x) + 1) +
  ylim(min(data_points$y) - 1, max(data_points$y) + 1) +
  guides(size = FALSE)
  

# Panel b ########################################################################################
# plotting dynamics
data_dynamics = data.frame(t = 0, biomass = 0, patch = 0)
for(k in 1:n_patch){
  temp_data = data.frame(t = out[,1], biomass = out[,k+1], patch = k)
  data_dynamics = rbind(data_dynamics, temp_data)
}
data_dynamics = data_dynamics[-1,]

panel_b = ggplot(data = data_dynamics) +
  geom_line(mapping = aes(x = t, y = biomass, group = patch)) +
  xlim(0, 60) +
  labs(x = "Time", y = "Biomass")+
  theme_classic()


# Panel c ########################################################################################
# plotting final state
panel_c = ggplot() +
  geom_line(data = data_links, mapping = aes(x = x, y = y, group = link, size = weight), color = "grey") +
  geom_point(data = data_points, mapping = aes(x = x, y = y, color = x_eq), size = 5) +
  labs(x="",y="",color="Biomass", size="Weight") +
  scale_size_continuous(range = c(0.05,1)) +
  scale_colour_gradient(low = "red", high ="green", limits = c(0, K*1.2)) +
  theme_classic() +
  labs(x="",y="") +
  theme(axis.text.x = element_blank(),legend.position="bottom",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "white")) +
  xlim(min(data_points$x) - 1, max(data_points$x) + 1) +
  ylim(min(data_points$y) - 1, max(data_points$y) + 1) +
  guides(size = FALSE)



ggarrange(panel_a, panel_b, panel_c, nrow = 1)



# one dimensional spectral reduction ##################################
odsr.object = odsr(W) # computing the one dimensional spectral reduction

odsr.object$a = sapply(odsr.object$a, Re)
odsr.object$alpha = Re(odsr.object$alpha)
odsr.object$alpha.out = Re(odsr.object$alpha.out)
odsr.object$beta = Re(odsr.object$beta)


# running the dimension reduction simulation:
odsr.out = odsr.dynamics(times = t,
                         x = x_init,
                         growth_function = logistic_allee,
                         dispersal_function = linear_dispersal,
                         odsr.object = odsr.object,
                         par = par)

#storing the output in a dataframe
data_odsr = data.frame(t = odsr.out[,1], R = odsr.out[,2])



## panel d #####################################################
dummy_data = data.frame(x = 0, y = 0, x_init = data_odsr$R[1], x_eq = data_odsr$R[length(t)])
panel_d = ggplot() +
  geom_point(data = dummy_data, mapping = aes(x = x, y = y, color = x_init), size = 25) +
  labs(x="",y="",color="") +
  scale_size_continuous(range = c(0.05,1)) +
  scale_colour_gradient(low = "red", high ="green", limits = c(0, K*1.2)) +
  theme_classic() +
  labs(x="",y="") +
  theme(axis.text.x = element_blank(),
        legend.position="",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "white")) +
  guides(size = FALSE)

## panel e ###########################################################
panel_e = ggplot() +
  geom_line(data = data_dynamics, mapping = aes(x = t, y = biomass, group = patch), alpha = 0.1) +
  geom_line(data = data_odsr, mapping = aes(x = t, y = R)) +
  xlim(0, 60) +
  labs(x = "Time", y = "Biomass")+
  theme_classic()

## panel_f ###########################################################

# deprecated
panel_f = ggplot() +
  geom_point(data = dummy_data, mapping = aes(x = x, y = y, color = x_eq), size = 25) +
  labs(x="",y="",color="") +
  scale_size_continuous(range = c(0.05,1)) +
  scale_colour_gradient(low = "red", high ="green", limits = c(0, K*1.2)) +
  theme_classic() +
  labs(x="",y="") +
  theme(axis.text.x = element_blank(),
        legend.position="",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "white")) +
  guides(size = FALSE)


panel_f_x = seq(0, 11, length.out = 1000)
panel_f_y = odsr.derivative(panel_f_x, panel_f_x, logistic_allee, dispersal_function = linear_dispersal, odsr.object, par)[[1]]
data_f = data.frame(x = panel_f_x, y = panel_f_y)


shift_axis <- function(p, y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
  
}

p = ggplot(data = data_f) +
  geom_point(mapping = aes(x = x, y = y), size = 0.1) +
  theme_classic() +
  labs(x = "R", y = "dR/dt") +
  scale_x_continuous(expand = c(0, 0), breaks = c(3, 6, 9), limits = c(0, 12))

panel_f = shift_axis(p, 0) + theme(axis.line.x = element_blank(), axis.title.x = element_text(vjust = 40, hjust = 0.98))


fig_1 = ggarrange(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, nrow = 2, ncol = 3, labels = letters[1:6])

ggsave("output/Figure_1.pdf",
       fig_1,
       height = 297/2,
       width = 210,
       units = "mm")

fig_1



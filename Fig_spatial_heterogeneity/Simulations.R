# a script that tests the one-dimensional spectral reduction over a list of networks

# sources and libraries ####################################################################
rm(list = ls())
library(deSolve)
library(igraph)
library(parallel)

# functions for the reduction
source("../one_dimension_spectral_reduction.R")

# functions for dynamics
source("../Dynamics/metapop_dynamics.R")
source("../Dynamics/Growth_functions.R")
source("../Dynamics/Dispersal_functions.R")
source("../Networks/Network_functions.R")


# creating output directories ###############################################################
dir.create("output")

# declaring parameters ######################################################################

# demographic properties
r = 1
baseline_K = 30
c = 0 # allee effect threshold == 0 --> the model is just a logistic growth


# main function, to loop over networks ######################################################

benchmark = function(par){
  network_id = par[1]
  k_sd = as.numeric(par[2])
  
  print(network_id)
  
  t = seq(0, 500, 1) # time seq for integration
  
  network = read_network(paste0("Network_list/", network_id))
  disp_matrix = network$dispersal_matrix
  n_patch = ncol(disp_matrix)
  
  
  data = data.frame(network_id = 0,
                    kind = network$summary$kind,
                    n_patch = 0,
                    rep = 0,
                    mu = 0,
                    k_sd = 0,
                    alpha = 0, alpha.out = 0, beta = 0, # odsr parameters
                    eigen.ratio = 0, # ratio of the first two eigen values
                    R_init = 0, # reduction of the initial state
                    R = 0, # predicted reduction
                    x_bar = 0, # average equilibrium
                    x_reduced = 0) # reduced eq
  
  for (rep in 1:3){ # we replicate different initial conditions for each network
    
    x_init = rnorm(n_patch, baseline_K, 0.1*baseline_K)
    K_vect = rnorm(n_patch, baseline_K, baseline_K * k_sd)
    K_vect = K_vect * (K_vect > 0)
    
    for (mu in seq(0.01, 1, length.out = 10)){# looping over mu
      W = disp_matrix*mu # scaling the disp_matrix by mu
      
      odsr.object = odsr(W) # computing the one dimensional spectral reduction
      
      if (Re(odsr.object$alpha) < 0 ){
        odsr.object = odsr(W, n_eig = 2)
      }
      
      odsr.object$a = sapply(odsr.object$a, Re)
      odsr.object$alpha = Re(odsr.object$alpha)
      odsr.object$alpha.out = Re(odsr.object$alpha.out)
      odsr.object$beta = Re(odsr.object$beta)
      
      par = list(r = 1,
                 K = K_vect,
                 c = 0)
      
      out = ode(y = x_init,
                times = t,
                func = metapop_dynamics,
                parms = par,
                growth_function = logistic_allee,
                dispersal_function = linear_dispersal,
                W = W)
      
      odsr.out = odsr.dynamics(times = t,
                               x = x_init,
                               growth_function = logistic_allee,
                               dispersal_function = linear_dispersal,
                               odsr.object = odsr.object,
                               par = par)
      
      R = odsr.out[nrow(odsr.out),2]
      
      eq_density = out[nrow(out), 2:(n_patch+1)]
      
      x_bar = mean(eq_density)
      x_reduced = t(odsr.object$a)%*%eq_density
      
      
      temp_data = data.frame(network_id = network_id,
                             kind = network$summary$kind,
                             n_patch = n_patch,
                             rep = rep,
                             mu = mu,
                             k_sd = k_sd,
                             alpha = odsr.object$alpha, alpha.out = odsr.object$alpha.out, beta = odsr.object$beta, # odsr parameters
                             eigen.ratio = abs(odsr.object$alpha/odsr.object$eigenvalues[2]), # ratio of the first two eigen values
                             R_init = t(odsr.object$a)%*%x_init, # reduction of the initial state
                             R = R, # predicted reduction
                             x_bar = x_bar, # average equilibrium
                             x_reduced = x_reduced) # reduced eq
      data = rbind(data, temp_data)
    }
  }
  write.table(data[-1,], paste("./output/", network_id, "_ksd", k_sd, ".csv", sep = ""))
  return(data[-1,])
}

# looping over all landscapes ######################################################################################################################

 ## generating parameters list ###########################################################################
  par = list()
  
  index = 1
  for(network_id in list.dirs('./Network_list', full.names = FALSE)[-1]){
    for(k_sd in seq(0, 0.15, length.out = 3)){
      par[[index]] = c(network_id, k_sd)
      index = index+1 
    }
  }

  ## looping in parallel ##################################################################################
  
  out = mclapply(par, benchmark, mc.cores = 6, mc.silent = F)
  
# Saving the output into a single .csv ####################################################################
  
  data = out[[1]]
  for (k in 2:length(out)){
    data = rbind(data, out[[k]])
  }
  write.csv(data, './output/0_All_data.csv')
  
  
  

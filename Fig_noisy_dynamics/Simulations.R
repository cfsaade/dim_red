# a script that tests the one-dimensional spectral reduction over a list of networks

# sources and libraries ####################################################################
rm(list = ls())
library(deSolve)
library(igraph)
library(parallel)
library(stats)

# functions for the reduction
source("../one_dimension_spectral_reduction.R")

# functions for dynamics
source("../Dynamics/metapop_dynamics.R")
source("../Dynamics/Growth_functions.R")
source("../Dynamics/Dispersal_functions.R")
source("../Networks/Network_functions.R")


# creating output directories ###############################################################
dir.create("output")

# declaring fixed parameters ######################################################################

# demographic properties
r = 1
baseline_K = 30
c = 0 # allee effect threshold == 0 --> the model is just a logistic growth

k_sd = 0
mu = 0.1


# timestep for the stochastic integration
dt = 0.01

#replicates:
n_rep = 3




# main function, to loop over networks ######################################################

benchmark = function(par){
  network_id = par[1]
  
  print(network_id)
  
  t = seq(0, 1000, 1) # time seq for integration
  
  network = read_network(paste0("Network_list/", network_id))
  disp_matrix = network$dispersal_matrix
  n_patch = ncol(disp_matrix)
  
  
  data = data.frame(network_id = '',
                    kind = network$summary$kind,
                    rep = 0,
                    
                    noise_amplitude = 0,
                    noise_cov_dist = 0,
                    
                    # odsr parameters
                    alpha = 0,
                    alpha.out = 0,
                    beta = 0,
                    odsr.cov = 0,
                    eigen.ratio = 0, # ratio of the first two eigen values
                    
                    R_init = 0, # reduction of the initial state
                    
                    # averages of the time series
                    R_mean = 0, # predicted reduction
                    x_bar_mean = 0, # average equilibrium
                    x_reduced_mean = 0, # reduced eq
                    
                    # sd-dev of the time series
                    R_sd = 0,
                    x_bar_sd = 0,
                    x_reduced_sd = 0,
                    
                    # acf lag 1 of time series
                    R_acf_1 = 0,
                    x_bar_acf_1 = 0,
                    x_reduced_acf_1 = 0)
  
  #for (rep in 1:3){ # we replicate different initial conditions for each network
   for (rep in 1:n_rep){ 
    x_init = rnorm(n_patch, baseline_K, 0.1*baseline_K)
    K_vect = rnorm(n_patch, baseline_K, baseline_K * k_sd)
    K_vect = K_vect * (K_vect > 0)
    
    for (noise_amplitude in seq(0.01, 0.1, length.out = 5)){ # looping over noise intensity
      for (dist_cov in c(0.05, 0.2, 0.4)){# looping over covariance distance
        network$cov_matrix = make_cov_matrix(network, dist_cov)*noise_amplitude # making a cov matrix of correlation distance dist_cov
        network$summary$covariance_distance = dist_cov
        W = disp_matrix*mu # scaling the disp_matrix by mu
        
        odsr.object = odsr(W, type = 'stochastic', cov_matrix = network$cov_matrix) # computing the one dimensional spectral reduction
        
        if (Re(odsr.object$alpha) < 0 ){
          odsr.object = odsr(W, type = 'stochastic', cov_matrix = network$cov_matrix, n_eig = 2)
        }
        
        odsr.object$a = sapply(odsr.object$a, Re)
        odsr.object$alpha = Re(odsr.object$alpha)
        odsr.object$alpha.out = Re(odsr.object$alpha.out)
        odsr.object$beta = Re(odsr.object$beta)
        odsr.object$cov = Re(odsr.object$cov)
        
        par = list(r = 1,
                   K = K_vect,
                   c = 0)
        
        out = ode(y = x_init,
                  times = t,
                  func = metapop_dynamics_stochastic,
                  parms = par,
                  growth_function = logistic_allee,
                  dispersal_function = linear_dispersal,
                  W = W,
                  
                  Cov = network$cov_matrix/dt,
                  method = "euler",
                  hini = dt)
        
        odsr.out = odsr.dynamics(times = t,
                                 x = x_init,
                                 growth_function = logistic_allee,
                                 dispersal_function = linear_dispersal,
                                 odsr.object = odsr.object,
                                 par = par,
                                 dt = dt)
        
        ## computing metrics for the 500 last timesteps of the dynamics:
          # extracting the last 500 timesteps
          x = out[500:nrow(out),2:ncol(out)]
          R = odsr.out[500:nrow(odsr.out), 2]
          
          #making x_bar and x_reduced time series:
          x_bar = rowMeans(x)
          x_reduced = x %*% odsr.object$a
          
          # averages of the time series
          R_mean = mean(R) # predicted reduction
          x_bar_mean = mean(x_bar) # average equilibrium
          x_reduced_mean = mean(x_reduced) # reduced eq
          
          # sd-dev of the time series
          R_sd = sd(R)
          x_bar_sd = sd(x_bar)
          x_reduced_sd = sd(x_reduced)
          
          # acf lag 1 of time series
          R_acf_1 = acf(R, 1, plot = F)[1]$acf[1]
          x_bar_acf_1 = acf(x_bar, 1, plot = F)[1]$acf[1]
          x_reduced_acf_1 = acf(x_reduced, 1, plot = F)[1]$acf[1]
          
        
        
        temp_data = data.frame(network_id = network_id,
                               kind = network$summary$kind,
                               rep = rep,
                               
                               noise_amplitude = noise_amplitude,
                               noise_cov_dist = dist_cov,
                               
                               # odsr parameters
                               alpha = odsr.object$alpha,
                               alpha.out = odsr.object$alpha.out,
                               beta = odsr.object$beta,
                               odsr.cov = odsr.object$cov,
                               eigen.ratio = abs(odsr.object$alpha/odsr.object$eigenvalues[2]), # ratio of the first two eigen values
                               
                               R_init = t(odsr.object$a)%*%x_init, # reduction of the initial state

                               # averages of the time series
                               R_mean = R_mean, # predicted reduction
                               x_bar_mean = x_bar_mean, # average equilibrium
                               x_reduced_mean = x_reduced_mean, # reduced eq
                               
                               # sd-dev of the time series
                               R_sd = R_sd,
                               x_bar_sd = x_bar_sd,
                               x_reduced_sd = x_reduced_sd,
                               
                               # acf lag 1 of time series
                               R_acf_1 = R_acf_1,
                               x_bar_acf_1 = x_bar_acf_1,
                               x_reduced_acf_1 = x_reduced_acf_1)
                               
                               
        data = rbind(data, temp_data)
      
      } # closing dist_cov loop
    } # closing noise_amplitude loop
  } # closing rep loop
  write.table(data[-1,], paste("./output/", network_id, "_ksd", k_sd, ".csv", sep = ""))
  return(data[-1,])
}

# looping over all landscapes ######################################################################################################################

 ## generating parameters list ###########################################################################
  par = list()
  
  index = 1
  for(network_id in list.dirs('./Network_list', full.names = FALSE)[-1]){
      par[[index]] = c(network_id)
      index = index+1 
  }

  ## looping in parallel ##################################################################################
  
  out = mclapply(par, benchmark, mc.cores = 6, mc.silent = F)
  
# Saving the output into a single .csv ####################################################################
  
  data = out[[1]]
  for (k in 2:length(out)){
    data = rbind(data, out[[k]])
  }
  write.csv(data, './output/0_All_data.csv')
  
  
  

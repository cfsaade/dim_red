# One dimension spectral reduction (ODSR) 

# Reduction function ######################################################################################

# function to call in order to get the parameters for an ODSR based on spectral analysis of the dispersal matrix
# (adapted from Laurence et Al.)

odsr = function(disp_matrix, type = 'deterministic', cov_matrix = 0, n_eig = 1){
  # odsr = One Dimensional Spectral Reduction
  # computes the quantities necessary for a one-dimensional reduction of a metapopulation
  #     based on spectral analysis of the matrix (adapted from Laurence et Al., 2019)
  #
  # Arguments:
  #   disp_matrix: dispersal matrix of the system
  #   type: whether to perform the odsr of a deterministic or stochastic system. values in ['deterministic', 'stochastic']
  #   cov_matrix: matrix of variance-covarianc of the noise for stochastic systems. 0 by default (dummy variable for 'deterministic' cases)
  #   n_eig: rank of the eigen-value for which to conduct the reduction (default = 1, for the most precise 1D reduction)
  
  if (type == 'deterministic'){
    out = odsr.deterministic(disp_matrix, n_eig)
    return(out)
  }
  else if (type == 'stochastic'){
    out = odsr.stochastic(disp_matrix, cov_matrix, n_eig)
    return(out)
  }
  else {
    stop(paste0("Unrecognised 'type' argument: '", type, "' instead of 'deterministic' or 'stochastic'."))
  }
}

# dynamics function ###########################################################

odsr.dynamics = function(times, x, growth_function, dispersal_function, odsr.object, par, dt = 0.001){
  # returns the dynamics of a system's 1-dimensional reduction.
  # Arguments:
  #   times: list of times at which the dynamics should be evaluated
  #   x: initial density in each patch
  #   growth_function: a function describing the local growth in each patch
  #   dispersal_function: a function describing how dispersal varies with species density
  #   odsr.object: the output of 'odsr' function when applied to the systems matrices (disp_ and cov_)
  #   par: a list of parameters to feed the growht and dispersal functions
  #   dt: integration time step, only useful for stochastic systems (deterministic systems are integrated using adaptive time steps)
  
  if (odsr.object$type == 'deterministic'){
    out = odsr.dynamics.deterministic(times, x, growth_function, dispersal_function, odsr.object, par)
  }
  else if (odsr.object$type == 'stochastic'){
    out = odsr.dynamics.stochastic(times, x, growth_function, dispersal_function, odsr.object, par, dt)
  }
  else {
    stop('Unrecognized odsr.object')
  }
  return(out)
}

# internal odsr functions #####################################################

  # functions called by 'odsr'
  
  ## deterministic odsr function ####################################################################################  

  odsr.deterministic = function(disp_matrix, n_eig = 1){
    # odsr = One Dimensional Spectral Reduction
    # computes the quantities necessary for a one-dimensional reduction of a metapopulation
    #     based on spectral analysis of the matrix (adapted from Laurence et Al., 2019)
    # disp_matrix is the dispersal matrix
    # n_eig is the rank of the eigen-value for which to conduct the reduction (default = 1, for the most precise 1D reduction)
    
    out = list()
    eig = eigen(t(disp_matrix))
    out$eigenvalues = eig$values
    
    #eigen value
    out$alpha = eig$values[n_eig]
    
    #normalized eigenvector
    out$a = eig$vectors[,n_eig]
    out$a = out$a/sum(out$a)
    
    #computing the vector of outgoing degrees:
    kout = apply(disp_matrix, 2, sum)
    
    #alpha.out = a^T times the vector of outgoing degrees
    out$alpha.out = t(out$a) %*% apply(disp_matrix, 2, sum)
    
    #beta
    out$beta = t(out$a) %*% (kout*diag(nrow(disp_matrix))) %*% out$a/(out$alpha.out * t(out$a) %*% out$a)
    
    
    out$type = 'deterministic'
    
    return(out)
  }
  
  ## Stochastic odsr function #################################################################################
  odsr.stochastic = function(disp_matrix, cov_matrix = 0, n_eig = 1){
    # odsr = One Dimensional Spectral Reduction
    # computes the quantities necessary for a one-dimensional reduction of a metapopulation
    #     based on spectral analysis of the matrix (adapted from Laurence et Al., 2019)
    # disp_matrix is the adjacency matrix
    # n_eig is the rank of the eigen-value for which to conduct the reduction (default = 1, for the most precise 1D reduction)
    # Cov: variance/covariance matrix in the case of noisy dynamics
    
    out = list()
    eig = eigen(t(disp_matrix))
    out$eigenvalues = eig$values
    
    #eigen value
    out$alpha = eig$values[n_eig]
    
    #normalized eigenvector
    out$a = eig$vectors[,n_eig]
    out$a = out$a/sum(out$a)
    
    #computing the vector of outgoing degrees:
    kout = apply(disp_matrix, 2, sum)
    
    #alpha.out = a^T times the vector of outgoing degrees
    out$alpha.out = t(out$a) %*% apply(disp_matrix, 2, sum)
    
    #beta
    out$beta = t(out$a) %*% (kout*diag(nrow(disp_matrix))) %*% out$a/(out$alpha.out * t(out$a) %*% out$a)
    
    #variance of resulting noise:
    out$cov = t(out$a) %*% cov_matrix %*% out$a
    
    out$type = 'stochastic'
    
    return(out)
  }
  
# internal dynamics functions ###################################################################
  
  ## deterministic dynamics functions ###########################################################
  
  # derivative
  odsr.derivative.deterministic = function(t, R, growth_function, dispersal_function, odsr.object, par){
    # computes the derivative of a single time step
    dR = growth_function(t, R, par)[[1]] + odsr.object$alpha*dispersal_function(t, R, par)[[1]] - odsr.object$alpha.out*dispersal_function(t, odsr.object$beta*R, par)[[1]]
    return(list(c(dR)))
  }
  

  
  odsr.dynamics.deterministic = function(times, x, growth_function, dispersal_function, odsr.object, par){
    # returns the dynamics of a deterministic system's 1-dimensional reduction.
    # Arguments:
    #   times: list of times at which the dynamics should be evaluated
    #   x: initial density in each patch
    #   growth_function: a function describing the local growth in each patch
    #   dispersal_function: a function describing how dispersal varies with species density
    #   odsr.object: the output of 'odsr' function when applied to the systems matrices (disp_ and cov_)
    #   par: a list of parameters to feed the growht and dispersal functions
    
    # reducing the initial density
    R = t(odsr.object$a)%*%x
    
    # reducing the parameters
    par.odsr = par
    for (k in 1:length(par)){
      par.odsr[[k]] = sum(odsr.object$a*par[[k]])
    }
    out = ode(y = R, times = times, func = odsr.derivative.deterministic,
              growth_function = growth_function,
              dispersal_function = dispersal_function,
              odsr.object = odsr.object,
              par = par.odsr)
    return(out)
  }
  
  ## stochastic dynamics functions ########################################################################################################################
  
  odsr.derivative.stochastic = function(t, R, growth_function, dispersal_function, odsr.object, par, dt){
    # computes the derivative of a single time step
    dR = growth_function(t, R, par)[[1]] + odsr.object$alpha*dispersal_function(t, R, par)[[1]] - 
      odsr.object$alpha.out*dispersal_function(t, odsr.object$beta*R, par)[[1]] +
      rnorm(1, 0, sqrt(odsr.object$cov/dt))*R
    return(list(c(dR)))
  }
  
  # odsr.object = odsr(disp_matrix)
  # odsr.derivative(0, 1, brummit, linear_dispersal, odsr.object, par = list(brummit.a = 0))
  
  odsr.dynamics.stochastic = function(times, x, growth_function, dispersal_function, odsr.object, par, dt){
    # returns the dynamics of a stochastic system's 1-dimensional reduction.
    # Arguments:
    #   times: list of times at which the dynamics should be evaluated
    #   x: initial density in each patch
    #   growth_function: a function describing the local growth in each patch
    #   dispersal_function: a function describing how dispersal varies with species density
    #   odsr.object: the output of 'odsr' function when applied to the systems matrices (disp_ and cov_)
    #   par: a list of parameters to feed the growht and dispersal functions
    
    # reducing the initial density
    R = t(odsr.object$a)%*%x
    
    # reducing the parameters
    par.odsr = par
    for (k in 1:length(par)){
      par.odsr[[k]] = sum(odsr.object$a*par[[k]])
    }
    out = ode(y = R, times = times, func = odsr.derivative.stochastic,
              growth_function = growth_function,
              dispersal_function = dispersal_function,
              odsr.object = odsr.object,
              par = par.odsr,
              dt = dt,
              
              method = "euler",
              hini = dt
    )
    return(out)
  }
  

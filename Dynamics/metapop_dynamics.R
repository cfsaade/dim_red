metapop_dynamics = function(t, x, growth_function, dispersal_function, W, par = c(), kout = apply(W, 2, sum)){
  # t = time (for the integrator)
  # x = vector of local densities
  # growth function: a shape for the local growth ("NoyMeir", "logistic_allee")
  # dispersal function: a shape for the dispersal ("linear_dispersal")
  # W: dispersal matrix
  # par: list of args to be passed to the growth and dispersal functions
  # kout = vector of outgoing links, computed from W by defautl
  
  local = growth_function(t, x, par)[[1]]
  spatial = dispersal_function(t, x, par)[[1]]
  dx = local - spatial*kout + W %*% spatial
  return(list(c(dx)))
}



metapop_dynamics_stochastic = function(t, x, growth_function, dispersal_function, W, par = c(), kout = apply(W, 2, sum), Cov){
  # t = time (for the integrator)
  # x = vector of local densities
  # growth function: a shape for the local growth ("NoyMeir", "logistic_allee")
  # dispersal function: a shape for the dispersal ("linear_dispersal")
  # W: dispersal matrix
  # par: list of args to be passed to the growth and dispersal functions
  # kout = vector of outgoing links, computed from W by default
  # Cov: var/cov matrix of the white noise. !!!! Should be scaled to the step size of the integrator !!!!
  
  
  local = growth_function(t, x, par)[[1]]
  spatial = dispersal_function(t, x, par)[[1]]
  white_noise = mvrnorm(1, rep(0, length(x)), Cov)*x
  
  
  dx = local - spatial*kout + W %*% spatial + white_noise
  return(list(c(dx)))
}


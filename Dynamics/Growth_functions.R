library(MASS)


brummit = function(t, x, par = list(brummit.a = 0, brummit.offset = 0, brummit.amplitude = 1)){
  ## a simple 3rd degree polynomial as used in Brummit et al. (2015) allowing for bistablility and catastrophic transitions 
  ## parameters:
  ##      a: control parameter, bistability for a in (xxx, xxx)
  ##      offset: parameter to shift the equilibriums to the left\right (the Brummit parameters result in a negative equilibrium)
  ##      amplitude: a parameter to make the derivative more (>1) or less (< 1) sharp, without changing the equilibria
  a = par$brummit.a
  offset = par$brummit.offset
  amplitude = par$brummit.amplitude
  
  x_eff = x-offset
  
  return(list(amplitude*(-x_eff**3 + x_eff + a)))
}


NoyMeir = function(t, n, par){
  # n = plant density
  r = par$r    # growth rate
  K = par$K    # carrying capacity
  B = par$B    # comsumption intensity
  A = par$A    # half consumption biomass
  
  dn = r*n*(1-n/K) - B*n/(A+n)
  
  return(list(dn))
}


logistic_allee = function(t, n, par){
  # logistic growth with an allee effect
  # n = species density density
  r = par$r    # growth rate
  K = par$K    # carrying capacity
  c = par$c    # critical density (in (0, 1), as a proportion of K): the density under which population growth is negative
  
  dn = r*n*(1-n/K)*(n/K - c)
  
  return(list(dn))
  
}


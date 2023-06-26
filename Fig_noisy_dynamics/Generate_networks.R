# generating networks


# sources and libraries ##############################################################
rm(list = ls())
library(igraph)

# functions for dynamics
source("../Networks/Network_functions.R")

# creating folders ###################################################################
dir.create("Network_list", recursive = T)


# Generating networks ################################################################

n_patch_list = c(21) # landscape sizes
n_landscapes = 10 # number of landscapes of each size
connectivity = 5 # connectivity for RGGs and erdos-renyi


  ## generating RGGs #################################################################
  id = 1
  for (n_patch in n_patch_list){
    for(k in 1:n_landscapes){
      print(paste0('RGG ', id))
      network = make_RGG(n_patch, connectivity)
      network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      while (is_connected(graph_from_adjacency_matrix(network$dispersal_matrix > 0)) != T){
        network = make_RGG(n_patch, connectivity)
        network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      }
      write_network(network, paste0('Network_list/RGG_', id))
      id = id+1
    }
  }
  
  ## generating OCNs #################################################################
  id = 1
  for (n_patch in n_patch_list){
    for(k in 1:n_landscapes){
      print(paste0('OCN  ', id))
      network = make_OCN(n_patch)
      network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      while (is_connected(graph_from_adjacency_matrix(network$dispersal_matrix > 0)) != T){
        network = make_OCN(n_patch)
        network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      }
      write_network(network, paste0('Network_list/OCN_', id))
      id = id+1
    }
  }

  ## generating erdos_renyi networks #################################################################
  id = 1
  for (n_patch in n_patch_list){
    for(k in 1:n_landscapes){
      print(paste0('ER  ', id))
      network = make_erdos_renyi(n_patch, c = connectivity)
      network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      while (is_connected(graph_from_adjacency_matrix(network$dispersal_matrix > 0)) != T){
        network = make_erdos_renyi(n_patch, c = connectivity)
        network$dispersal_matrix = randomize_link_strength(network$dispersal_matrix, 0.2)
      }
      write_network(network, paste0('Network_list/ER_', id))
      id = id+1
    }
  }


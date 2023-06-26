# generating networks


# sources and libraries ##############################################################
rm(list = ls())
library(igraph)

# functions for dynamics
source("../Networks/Network_functions.R")

# generating the networks ############################################################

OCN = make_OCN(21)
plot_network(OCN)
write_network(OCN, 'OCN')


RGG = make_RGG(21, 4)
plot_network(RGG)
write_network(RGG, 'RGG')


ER = make_erdos_renyi(21, 4)
plot_network(ER)
write_network(ER, 'ER')



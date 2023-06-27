Scripts and data for "Spectral dimension reduction of metapopulation dynamics"
Camille Saade, Emanuel A. Fronhofe, Sonia Kéfi

Abstract:

In nature, species interact with each other in spatially structured landscapes. 
Modelling such landscapes often relies on multidimensional dynamical systems, where one equation describes a species density at a given location. With such a formalism, the complexity of the model grows with the number of species and locations in the ecosystem described. The mathematical analysis of such systems becomes quickly difficult for more than a few species or habitable patches. Focusing on spatial complexity, we develop an approach allowing to simplify the dynamics of spatially structured populations. We consider a single species living in a patchy environment and show that its dynamics --- described by $N$ coupled differential equations each describing a patch --- can be reduced to a single one-dimensional equation. The reduction approach captures the overall state of the metapopulation, predicting its average biomass and temporal variability. We demonstrate its efficiency on different landscape structures and in various settings, including under spatially heterogeneous environmental conditions and with environmental temporal stochasticity. This approach greatly reduces the complexity of spatially structured systems, reducing their properties to a few parameters, and allows a better understanding of metapopulation dynamics.

main scripts and utility functions:
- one_dimension_spectral_reduction.R : script with various functions to conduct the dimension reduction of a given spatial landscape.
- Dynamics : folder with scripts related to the dynamics (dispersal functions, growth functions and metapop dynamics).
- Networks : utility functions for network generation and manipulations.
    
Folders named "Fig_xxx" are each responsible for running a set of simulations and plot the associated figure(s).
In each folder, run locally Generate_networks.R, then Simulations.R, then plot_figures.R to reconduct the full analysis.


Plot_networks is a folder to plot some example networks for illustration purposes.

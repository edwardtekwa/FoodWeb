# FoodWeb
Matlab foodweb simulation in a spatial temperature gradient with temperature change over time
Edward Tekwa Oct 26, 17

Important files:

Slurm_script2.m: run simulation script on multiple processors on cluster
Make_warming_MultPtsstats_Parallel0.m: run a set of simulations, with each running a burn-in period after which a no-warming scenario and a warming scenario diverge.  Only stats from the end of the simulations are saved in a .mat file.
make_parameters.m: specify how long to run each simulation, how many degrees warming
make_traits.m: specify biological parameters in arrays, in particular: sdm, sdv, td, basalSize, as well as general parameters such as number of species, number of patches, and environmental and biological temperature ranges
sub_demog.m: function called during simulation to determine who eats what
sub_feedfunc.m: function called from make_parameters.m to specify the vulnerability matrix
sub_move.m: function called during simulation to determine who moves where
Add_warming_endstats_movement_fixedDir_new.m: run after a set of simulations are saved as .mat to collect data and plot basic statistics
plot_demog_spatial_endSnap.m: first load a .mat file of a single simulation, then plot time series and spatial distributions

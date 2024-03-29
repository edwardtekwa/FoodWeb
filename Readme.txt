Readme for https://github.com/edwardtekwa/FoodWeb
Tekwa, Watson, Pinsky, Body-size and food-web interactions mediate species range shifts under warming

Overview: Code for simulating food web dynamics under warming and plotting results. Runs on Matlab versions R2019a to R2020a.

Main files (in order of appearance in Instructions):
skewThEnv.m
Slurm_script_Foodweb.sh
Make_warming_MultPtsstats_Parallel.m
make_traits.m
make_parameters.m
sub_feedfunc_adj.m
Add_warming_endstats_dispersal_fixedDir.m
estSingleSpeciesModel.m
Add_warming_endstats_dispersal_plots.m
WarmingDispersalStats_main.mat
plot_demog_spatial_endS.m
Foodweb_numSpecies200_dT3_basalSize0.01_meanD3_pInedible0_fIII.mat

Load all files into one directory. Auxiliary files are called by the main files.

Instructions:

1. To run simulations from scratch on a remote cluster, run Slurm_script_Foodweb.sh after editing information on cluster nodes. Otherwise, to run simulations from scratch on a local computer, run Make_warming_MultPtsstats_Parallel.m after editing parallel computing parameters at the top. Output are 10 replicates for each of 6 dispersal rates in .mat files. NOTE: it takes about 23 core days on a 2020 Intel Xeon processor to complete this set.

Edit Make_warming_MultPtsstats_Parallel.m to change the number of replicates per dispersal rate.
Edit make_traits.m to specify warming rate, initial conditions, and main spatial and biological parameters that can take multiple experimental values.
Edit make_parameters.m to specify biological parameters, i.e. for sensitivity tests.
Edit sub_feedfunc_adj.m to specify feeding parameters.

2. To process data from new simulations and run single-species projections, run Add_warming_endstats_dispersal_fixedDir.m after editing "Path" to the local folder containing the simulations.

Edit estSingleSpeciesModel.m to change parameters for the single-species counterfactual model fit to food web data in the no-warming period. Data for subsequent plotting is saved as WarmingDispersalStats_<TimeData>.mat.

3. To plot aggregate results, run Add_warming_endstats_dispersal_plots.m after editing the .mat file target under "load". Otherwise, run as is to plot data from the paper's main simulations (in WarmingDispersalStats_main.mat).

4. To plot time series and food web diagram for one simulation replicate, run plot_demog_spatial_endS.m. This plots data from a sample replicate (Foodweb_numSpecies200_dT3_basalSize0.01_meanD3_pInedible0_fIII.mat).

Edit the .mat file target under "load" to use another simulation replicate. Edit patches and time points for food web diagrams.

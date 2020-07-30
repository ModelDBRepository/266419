This code contains conductance-based models of Dopaminergic (DA) and GABAergic neurons, used in "Distinct temporal structure of nicotinic ACh receptor activation determines responses of VTA neurons to endogenous ACh and nicotine" eNeuro paper.

The model is presented in MATLAB format (mainNicandAChDAmodel.m) with associated built MEX files and source C code (Both_NicandAChDAmodel.cpp).

mainNicandAChDAmodel.m is the main program, which calculates the firing rate and %spikes within bursts of dopamine neuron in 5 cases: knock-out of β2-containing nAChRs, β2-containing nAChRs only on DA neurons, the nAChRs only on GABA neurons, the nAChRs on both DA and GABA neurons and “wild” type. This script reproduces figures 5 and 6 of the manuscript. Generated plots might look insignificantly different from the published figures as the spikes for ACh and Glu inputs are drawn from the bimodal and Poisson distributions respectively for each trial.

This code contains complied MEX files for 64-bit system and the source C code that could be compiled on your platform.  We recommend compiling MEX files on your platform. This can be done by uncommenting the "Both_NicandAChDAmodel.cpp " line (line 1) in the mainNicandAChDAmodel.m (requires MEX compiler installed within MATLAB).

Both_NicandAChDAmodel.cpp is a C script for calculating membrane potential of the DA neuron.

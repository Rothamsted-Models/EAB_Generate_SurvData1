The Emerald Ash Borer (EAB; Agrilus planipennis Fairmare) is a destructive pest for ash trees causing the death of millions of trees in parts of the globe where it is non-native. This pest poses a threat to ash in Great Britain (GB) and so establishing well-defined surveillance strategies to detect EAB is essential. We have developed a model of EAB lifecycle and spread across the Great Britain. The model outputs can then be used to predict the best places to locate surveillance technologies (e.g., girdled trees or traps).

The model is written in C++ and was compiled using Microsoft Visual Studio 2022. Stochastic sampling and integration of the dispersal kernel use NAG libraries https://nag.com/nag-library/  but those calls can be replaced by other libraries if you do not have a license. 

The model requires the entry point to be initialized. We developed maps to identify the most likely first incursion point at locations across Great Britain on a 1km x 1km grid. These maps are informed by data from the forestry commission on likely entry points and data on firewood use (see below). When creating these maps, we considered 9 scenarios related to how certain we are that EAB will arrive through known pathways related to wood imports (70%, 50%, 30%) and the probability that EAB would escape at port rather than at the onwards depots (25%, 50%, 75%).  For each scenario we sampled 10000 realizations for likely entry points.  These form inputs to our model “SampleEntry.txt”.  

The model is set up to run for annual time steps for 8 years. For each simulated year the larval density and tree health is recorded in infested cells and reported in “DataForOpt.txt”. The model is currently set up to do this for up to 10000 realisations of entry point. The data in the output file can then be used to optimise surveillance. 

This code was produced as part of the NERC funded Smarties project NE/T007729/1. Rothamsted Research receives strategic funding from the Biotechnology and Biological Sciences Research Council of the United Kingdom (BBSRC). We also acknowledge support from the Growing Health Institute Strategic Programme [BB/X010953/1; BBS/E/RH/230003C]

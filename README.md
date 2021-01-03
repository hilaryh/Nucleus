# Nucleus
Processing/analysing data

cutNuc1DData.m runs the scripts to load, smooth, and analyse data in order. Paths to read experimental data and metadata are set up for my computer and file types.
	homeDir/data_montage should point to where data files are stored
	processDir should point to where plots should be saved
	There may be references to a saveDir in some functions that I've missed. This can point to the same location as processDir to save all plots in the same folder.
	
Data is smoothed with butterworthImg.m and processImg.m
	processImg.m is tailored to a specific 1D calcium fluorescence dataset. If it does not appear to work on your dataset, you may need to experiment with the threshold variable 'minProm' on line 13 (which describes a minimum difference between the average calcium fluorescence and fluorescence at a peak for it to count as a calcium transient) and the inputs to the function 'flattenInit' on line 23 (In this line 'plocs(1)-40' is an time point estimated to occur before any calcium transients have been initiated)
	Additionally, maxTime is the last time point that this function will look at. Depending on the resolution and length of your dataset, this number may also need to be changed.
	
findScanBorder.m estimates the border of the cell in each linescan. 
	If this function does not appear to be working, try changing the variable 'timePoint' to a time point soon after a calcium transient for maximum contrast between the cell and background fluorescence.
	'cutOff', the threshold fluorescence above which the function assumes it is looking at a cell may also need to be adjusted. To visualise what this function is doing, run it with the second input variable (redoBool) equal to 1.

findCriticalPoints.m estimates the location of each calcium transient. 
	Within the function, 'peakSeparation' is the minimum time between transients. If your transients are occuring closer together, or much further apart, than 1500 pixels, change this value to something more approptiate.
	
findNuc.m attempts to identify the nucleus in each dataset based on the difference in timing between calcium transient initiation and the subsequent peak at each point along the linescan. This can be fairly finicky. I found that it was generally only accurate to within a few pixels so, depending on the size of your dataset, it may be easier to hard code the location of the nucleus in each cell.
	
Splines are approximated from data with avCyt_spline.m 

nonlinearDiff.m solves a numerical, 1D model of Ca2+ diffusion into the nucleus. It takes as input the nuclear diffusion coefficient, an array containing buffer and nuclear size parameters, parameters for a spline approximation of the cytosolic Ca2+ transient, interpolation points of the spline approximation, and the temporal resolution, and simulation time.

sim_nucleus_cyl_spline.m solves the analytical solution to the same model (numerically, with infinite sums truncated at maxt=300).

sim_nucleus_cyl_wz_spline_no_cp_c.m solves an analytical 3D model of Ca2+ diffusion into the nucleus.

To generate plots exploring the numerical model, run plotNumDiffPS.m
To generate plots exploring and comparing the analytical models, run plotAnaNuc.m

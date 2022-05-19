# sgd_shuffling
Experiments for the paper "Computing the Variance of Shuffling Stochastic Gradient Algorithms via Power Spectral Density Analysis"

Given configurations of n (number of functions), d (dimension), gamma (stepsize), alpha (momentum parameter) and runs 
(number of runs to get average and standard deviations), run 'sgd_averages.m', 'sgdm_averages.m' and 'snag_averages.m'.
'sgd_averages.m' runs SGD with replacement, SGD-RR and SGD-SO, each under the standard stochastic noise and under the
zero-th order noise model introduced in the paper. For each algorithm, it prints the mean squared error averaged over runs
and the standard deviation of the estimate, which have been used to populate the tables of the paper. It also keeps in the
the MATLAB workspace the sequence of distances to the optimum for the last run of each algorithm, under variable names like
iterates_dist_SGD_RR or iterates_dist_SGD_SO_approx. 'sgdm_averages.m' and 'snag_averages.m' are analogous for SGDM and 
SNAG, respectively.

Once these scripts have been run and the variables iterates_dist_* are stored in the workspace, the figures can be produced 
by executing 'plots_paper_tiled.m', 'plots_paper_tiled_2.m' and 'plots_paper_tiled_3.m', which were used to obtain Figures 3,
1 and 2 of the paper, respectively.

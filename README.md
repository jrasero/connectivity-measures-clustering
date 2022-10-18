# connectivity-measures-clustering

Repository that clusters two different functional connectivity measures and four time series extraction strategies.

The code is as follows:

- "compute_similarities_whole_brain.m" : For each fMRI session and frequency band, it first concatenates the connectivity maps across seed regions 
and then uses this big vector to compute the 8 x 8 similarity matrix  (2 connectivity measures, 4 TS extraction strategies).

- "compute_similarities_per_roi.m" : Same as above, but here a similarity matrix is computed for each seed region, using their hemisphere-wise concatenated maps.

- "perform_clustering.m" : For each band separately, we first average similarity matrices across sessions. Then, a distance matrix is computed as 1-R, where
R is the average similarity matrix in each scenario (i.e. whole-brain or individual region). Finally, a hierarchical clustering with an 'average' linkage method
is performed to see how functional connectivity measures and time series extraction pipelines group together.



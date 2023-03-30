# connectivity-measures-clustering

Repository that clusters two different functional connectivity measures (PLV and ciPLV) and four time series extraction strategies (MaxAfter, MeanAfter, MeanBefore and PCABefore).

The code is organised as follows:

- "compute_similarities.m": This is just the script that runs the following three scripts in the given order to compute the similarity matrices between pipelines:

  1. "compute_similarities_whole_brain.m" : For each fMRI session and frequency band, it first concatenates the connectivity maps across seed regions 
  and then uses this big vector to compute the 8 x 8 similarity matrix  (2 connectivity measures, 4 TS extraction strategies).

  2. "compute_similarities_per_roi.m" : Same as above, but here a similarity matrix is computed for each seed region, using their hemisphere-wise concatenated maps.

  3. "compute_similarities_per_roi_and_hemisphere.m" : Same as above, but here a similarity matrix is computed for each seed region and hemisphere sepaarately"

- "perform_clustering.m" : It computes the clustering as follows. for each band separately, it first averages the similarity matrices across sessions. Then, a distance matrix is computed as 1-R, where  R is the average similarity matrix in each scenario (i.e. whole-brain or individual region, see below). Finally, a hierarchical clustering with an 'average' linkage method is performed to see how functional connectivity measures and time series extraction pipelines group together. This script runs the following three scripts in the given order to perform the clustering between pipelines:

  1. "clustering_whole_brain.m" : It performs the clustering for the whole-brain.
  2. "clustering_roi_34.m" : It performs the clustering for the regions with left and right hemispheres together.
  3. "clustering_roi_68.m" : It performs the clustering for the regions and left and right hemispheres separately.
  
  



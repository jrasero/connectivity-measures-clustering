clear; clc;

% Whole brain case
run("compute_similarities_whole_brain.m")
% 34 regions (left and right hemispheres merged) case
run("compute_similarities_per_roi.m")
% 68 regions (hemispheres separately) case
run("compute_similarities_per_roi_and_hemisphere.m")
  

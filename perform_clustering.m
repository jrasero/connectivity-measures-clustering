clear; clc;

% Whole-brain case
run("clustering_whole_brain.m")
% 34 regions (left and right hemispheres merged) case
run("clustering_roi_34.m")
% 68 regions (hemispheres separately) case
run("clustering_roi_68.m")
  

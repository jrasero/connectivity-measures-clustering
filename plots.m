clc;clear;
%% Correlate distances with regions sizes

% load recording names from previously save mat file
load('pipeline_options.mat')

n_measures = length(measures);
n_funcs = length(funcs);
n_bands = length(bands);
n_rois = length(rois_merged);
n_recs = 160;

lh_vols = readtable("data/cvs_avg35_inMNI152_desikan_vols_lh.txt");
rh_vols = readtable("data/cvs_avg35_inMNI152_desikan_vols_rh.txt");

% Create volumes (34 regions, merged hemispheres)
vols_34 = {};
for iroi=1:n_rois
    tic
    roi = rois_merged{iroi};
    disp(['Doing ROI:' ' ' roi])
    vols_34.(roi) = lh_vols.(strcat('lh_', roi, '_volume')) + rh_vols.(strcat('rh_', roi, '_volume'));    
end

% Create volumes (68 regions, separate hemispheres)
vols_68 = {};
for iroi=1:n_rois
    tic
    roi = rois_merged{iroi};
    disp(['Doing ROI:' ' ' roi])
    vols_68.(strcat(roi, "_L")) = lh_vols.(strcat('lh_', roi, '_volume'));
    vols_68.(strcat(roi, "_R")) = rh_vols.(strcat('rh_', roi, '_volume'));
end
        
band = 'alpha';

d_alpha_34 = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('clustering/roi_34/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_alpha_34.(roi) = mean(squareform(D));
end

d_alpha_68 = {};
fieldNames = fieldnames(vols_68);
for iroi=1:numel(fieldNames)
    roi = fieldNames{iroi};
    output_dir = strcat('clustering/roi_68/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_alpha_68.(roi) = mean(squareform(D));
end

band = 'beta';

d_beta_34 = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('clustering/roi_34/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_beta_34.(roi) = mean(squareform(D));
end

d_beta_68 = {};
fieldNames = fieldnames(vols_68);
for iroi=1:numel(fieldNames)
    roi = fieldNames{iroi};
    output_dir = strcat('clustering/roi_68/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_beta_68.(roi) = mean(squareform(D));
end

band = 'theta';

d_theta_34 = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('clustering/roi_34/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_theta_34.(roi) = mean(squareform(D));
end

d_theta_68 = {};
fieldNames = fieldnames(vols_68);
for iroi=1:numel(fieldNames)
    roi = fieldNames{iroi};
    output_dir = strcat('clustering/roi_68/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_theta_68.(roi) = mean(squareform(D));
end

% Plot distances across bands
f = figure('visible','off');
boxplot(transpose([struct2array(d_alpha_34); struct2array(d_beta_34); struct2array(d_theta_34)]),...
    'Notch','on','Labels', {'Alpha','Beta', 'Theta'})
title('Distances distribution')
saveas(f, 'plots/dist_vs_bands_34.png')
close(f)

% Plot distances across bands
f = figure('visible','off');
boxplot(transpose([struct2array(d_alpha_68); struct2array(d_beta_68); struct2array(d_theta_68)]),...
    'Notch','on','Labels', {'Alpha','Beta', 'Theta'})
title('Distances distribution')
saveas(f, 'plots/dist_vs_bands_68.png')
close(f)

% Compute correlation between distances and Volumes
% Alpha
[R1, P1] = corrcoef(struct2array(vols_34), struct2array(d_alpha_34));

f = figure('visible','off');
scatter(struct2array(vols_34), struct2array(d_alpha_34), 'b','*')
title(strcat('Alpha Band: R = ', num2str(R1(1,2)), ',', ' P = ', num2str(P1(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_alpha_34.png')
close(f)

[R1, P1] = corrcoef(struct2array(vols_68), struct2array(d_alpha_68));

f = figure('visible','off');
scatter(struct2array(vols_68), struct2array(d_alpha_68), 'b','*')
title(strcat('Alpha Band: R = ', num2str(R1(1,2)), ',', ' P = ', num2str(P1(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_alpha_68.png')
close(f)

% Beta
[R2, P2] = corrcoef(struct2array(vols_34), struct2array(d_beta_34));

f = figure('visible','off');
scatter(struct2array(vols_34), struct2array(d_beta_34), 'mo')
title(strcat('Beta Band: R = ', num2str(R2(1,2)), ',', ' P = ', num2str(P2(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_beta_34.png')
close(f)

[R2, P2] = corrcoef(struct2array(vols_68), struct2array(d_beta_68));

f = figure('visible','off');
scatter(struct2array(vols_68), struct2array(d_beta_68), 'mo')
title(strcat('Beta Band: R = ', num2str(R2(1,2)), ',', ' P = ', num2str(P2(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_beta_68.png')
close(f)

% Theta
[R3, P3] = corrcoef(struct2array(vols_34), struct2array(d_theta_34));

f = figure('visible','off');
scatter(struct2array(vols_34), struct2array(d_theta_34), 'rx')
title(strcat('Theta Band: R = ', num2str(R3(1,2)), ',', ' P = ', num2str(P3(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_theta_34.png')
close(f)

[R3, P3] = corrcoef(struct2array(vols_68), struct2array(d_theta_68));

f = figure('visible','off');
scatter(struct2array(vols_68), struct2array(d_theta_68), 'rx')
title(strcat('Theta Band: R = ', num2str(R3(1,2)), ',', ' P = ', num2str(P3(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_vols_theta_68.png')
close(f)
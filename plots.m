clc;clear;
%% Correlate distances with regions sizes

% load recording names from previously save mat file
load('pipeline_options.mat')

n_measures = length(measures);
n_funcs = length(funcs);
n_bands = length(bands);
n_rois = length(rois_merged);
n_recs = 160;

[num,txt,raw] = xlsread("data/Seed_Size.xlsx");

seed_size={};
for ii=1:68
    seed_size.(strcat(string(raw(ii,1)), "_", string(raw(ii,2)))) = num(ii);
end
    

% lh_vols = readtable("data/cvs_avg35_inMNI152_desikan_vols_lh.txt");
% rh_vols = readtable("data/cvs_avg35_inMNI152_desikan_vols_rh.txt");
% 
% % Create volumes (34 regions, merged hemispheres)
% vols_34 = {};
% for iroi=1:n_rois
%     tic
%     roi = rois_merged{iroi};
%     disp(['Doing ROI:' ' ' roi])
%     vols_34.(roi) = lh_vols.(strcat('lh_', roi, '_volume')) + rh_vols.(strcat('rh_', roi, '_volume'));    
% end
% 
% % Create volumes (68 regions, separate hemispheres)
% vols_68 = {};
% for iroi=1:n_rois
%     tic
%     roi = rois_merged{iroi};
%     disp(['Doing ROI:' ' ' roi])
%     vols_68.(strcat(roi, "_L")) = lh_vols.(strcat('lh_', roi, '_volume'));
%     vols_68.(strcat(roi, "_R")) = rh_vols.(strcat('rh_', roi, '_volume'));
% end

size_34 = {};
for iroi=1:n_rois
    tic
    roi = rois_merged{iroi};
    disp(['Doing ROI:' ' ' roi])
    size_34.(roi) = seed_size.(strcat(roi, '_L')) + seed_size.(strcat(roi, '_R'));    
end

% Create volumes (68 regions, separate hemispheres)
size_68 = {};
for iroi=1:n_rois
    tic
    roi = rois_merged{iroi};
    disp(['Doing ROI:' ' ' roi])
    size_68.(strcat(roi, "_L")) = seed_size.(strcat(roi, '_L'));
    size_68.(strcat(roi, "_R")) =  seed_size.(strcat(roi, '_R'));
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
fieldNames = fieldnames(size_68);
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
fieldNames = fieldnames(size_68);
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
fieldNames = fieldnames(size_68);
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
saveas(f, 'plots/dist_vs_bands_34_matlab.png')
close(f)

% Plot distances across bands
f = figure('visible','off');
boxplot(transpose([struct2array(d_alpha_68); struct2array(d_beta_68); struct2array(d_theta_68)]),...
    'Notch','on','Labels', {'Alpha','Beta', 'Theta'})
title('Distances distribution')
saveas(f, 'plots/dist_vs_bands_68_matlab.png')
close(f)

% Compute correlation between distances and seed size
% Alpha
[R1, P1] = corrcoef(struct2array(size_34), struct2array(d_alpha_34));

f = figure('visible','off');
scatter(struct2array(size_34), struct2array(d_alpha_34), 'b','*')
title(strcat('Alpha Band: R = ', num2str(R1(1,2)), ',', ' P = ', num2str(P1(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_alpha_34_matlab.png')
close(f)

[R1, P1] = corrcoef(struct2array(size_68), struct2array(d_alpha_68));

f = figure('visible','off');
scatter(struct2array(size_68), struct2array(d_alpha_68), 'b','*')
title(strcat('Alpha Band: R = ', num2str(R1(1,2)), ',', ' P = ', num2str(P1(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_alpha_68_matlab.png')
close(f)

% Beta
[R2, P2] = corrcoef(struct2array(size_34), struct2array(d_beta_34));

f = figure('visible','off');
scatter(struct2array(size_34), struct2array(d_beta_34), 'mo')
title(strcat('Beta Band: R = ', num2str(R2(1,2)), ',', ' P = ', num2str(P2(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_beta_34_matlab.png')
close(f)

[R2, P2] = corrcoef(struct2array(size_68), struct2array(d_beta_68));

f = figure('visible','off');
scatter(struct2array(size_68), struct2array(d_beta_68), 'mo')
title(strcat('Beta Band: R = ', num2str(R2(1,2)), ',', ' P = ', num2str(P2(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_beta_68_matlab.png')
close(f)

% Theta
[R3, P3] = corrcoef(struct2array(size_34), struct2array(d_theta_34));

f = figure('visible','off');
scatter(struct2array(size_34), struct2array(d_theta_34), 'rx')
title(strcat('Theta Band: R = ', num2str(R3(1,2)), ',', ' P = ', num2str(P3(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_theta_34_matlab.png')
close(f)

[R3, P3] = corrcoef(struct2array(size_68), struct2array(d_theta_68));

f = figure('visible','off');
scatter(struct2array(size_68), struct2array(d_theta_68), 'rx')
title(strcat('Theta Band: R = ', num2str(R3(1,2)), ',', ' P = ', num2str(P3(1,2))))
xlabel('Seed size (# Vertices)')
ylabel('Average distance')
lsline
saveas(f, 'plots/dist_vs_size_theta_68_matlab.png')
close(f)
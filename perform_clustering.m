% load recording names from previously save mat file
load('pipeline_options.mat')

n_measures = length(measures);
n_funcs = length(funcs);
n_bands = length(bands);
n_rois = length(rois_merged);
n_recs = 160;

ticklabels = {};
kk=1;
for ii=1:n_measures
    for jj=1:n_funcs
        ticklabels{kk} = strcat(measures{ii}, '\_', funcs{jj});
        kk = kk + 1;
    end
end

% Iterate through bands
 for iband=1:n_bands
    tic
    band = bands{iband};
    disp(['Doing band:' ' ' band])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% First whole-brain %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    mkdir(strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',... 
        'whole_brain'), band)

    output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',...
                   'whole_brain', "/", band);

    % 1- Load each similarity matrix per recording and compute its average
    R_rec = zeros(n_recs, n_measures*n_func, n_measures*n_func);
    for irec=1:n_recs
        load(strcat('similarity_matrices_roi/whole_brain', ...
            "/", band, "/", 'recording_', num2str(irec), '.mat'), 'R')

        R_rec(irec, :, :) = R;
    end
    R_avg = squeeze(mean(R_rec, 1)); % average across recordings

    % 2-Plot this matrix
    f = figure('visible','off');
    imagesc(R_avg)
    xticks(1:n_measures*n_func)
    xticklabels(ticklabels)
    xtickangle(45)
    yticks(1:n_measures*n_func)
    yticklabels(ticklabels)
    colorbar()
    saveas(f, strcat(output_dir, '/', 'average_similarity_matrix.png'))
    close(gcf)

    % 3-Apply hierarchical clustering directly to distance matrix from
    % from the similarity matrix

    D = 1-R_avg; 

    Z = linkage(squareform(D), 'average');
    f = figure('visible','off');
    dendrogram(Z, 'Labels', ticklabels)
    title("Hierarchical clustering", 'fontsize', 20)
    saveas(gcf, strcat(output_dir, '/', 'hierarchical_clustering.png'))
    close(gcf)
    save(strcat(output_dir, '/', 'hierarchical_clustering.mat'),...
        'Z', 'D', 'R_avg')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Each ROI separately %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iroi=1:n_rois
        tic
        roi = rois_merged{iroi};
        mkdir(strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',... 
            roi), band)
        output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',...
            roi, "/", band);
        R_rec = zeros(n_recs, n_measures*n_func, n_measures*n_func);
        
        for ii=1:n_recs
            load(strcat('similarity_matrices_roi', "/",...
                roi, "/", band, "/",  ...
                'recording_', num2str(irec), '.mat'), 'R')
              R_rec(ii, :, :) = R;
        end
        R_avg = squeeze(mean(R_rec, 1)); % average across recordings
          
         % 2-Plot this matrix
         f = figure('visible','off');
         imagesc(R_avg)
         xticks(1:n_measures*n_func)
         xticklabels(ticklabels)
         xtickangle(45)
         yticks(1:n_measures*n_func)
         yticklabels(ticklabels)
         colorbar()
         saveas(f, strcat(output_dir, '/', 'average_similarity_matrix.png'))
         close(f)
          
         % 3-Apply hierarchical clustering directly to distance matrix from
         % from the similarity matrix
         D = 1-R_avg; 
         Z = linkage(squareform(D), 'average');
         f = figure('visible','off');
         dendrogram(Z, 'Labels', ticklabels)
         title("Hierarchical clustering", 'fontsize', 20)
         saveas(gcf, strcat(output_dir, '/', 'hierarchical_clustering.png'))
         close(gcf)
         save(strcat(output_dir, '/', 'hierarchical_clustering.mat'),...
             'Z', 'D', 'R_avg')
    end
    toc
 end

%% Correlate distances with regions sizes
lh_vols = readtable("cvs_avg35_inMNI152_desikan_vols_lh.txt");
rh_vols = readtable("cvs_avg35_inMNI152_desikan_vols_rh.txt");

vols = {};
for iroi=1:n_rois
    tic
    roi = rois_merged{iroi};
    disp(['Doing ROI:' ' ' roi])
    vols.(roi) = lh_vols.(strcat('lh_', roi, '_volume')) + rh_vols.(strcat('rh_', roi, '_volume'));    
end
         
band = 'alpha';
d_alpha = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_alpha.(roi) = mean(squareform(D));
end

band = 'beta';
d_beta = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_beta.(roi) = mean(squareform(D));
end

band = 'theta';
d_theta = {};
for iroi=1:n_rois
    roi = rois_merged{iroi};
    output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/',...
        roi, "/", band);
    load(strcat(output_dir, '/', 'hierarchical_clustering.mat'), 'D')
    d_theta.(roi) = mean(squareform(D));
end

% Plot distances across bands
f = figure('visible','off');
boxplot(transpose([struct2array(d_alpha); struct2array(d_beta); struct2array(d_theta)]),...
    'Notch','on','Labels', {'Alpha','Beta', 'Theta'})
title('Distances distribution')
saveas(f, ...
    '/home/javi/Documentos/measures_clustering_daniele/clustering/dist_vs_bands.png')
close(f)

% Compute correlation between distances and Volumes
% Alpha
[R1, P1] = corrcoef(struct2array(vols), struct2array(d_alpha));

f = figure('visible','off');
scatter(struct2array(vols), struct2array(d_alpha), 'b','*')
title(strcat('Alpha Band: R = ', num2str(R1(1,2)), ',', ' P = ', num2str(P1(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, ...
    '/home/javi/Documentos/measures_clustering_daniele/clustering/dist_vs_vols_alpha.png')
close(f)

% Beta
[R2, P2] = corrcoef(struct2array(vols), struct2array(d_beta));

f = figure('visible','off');
scatter(struct2array(vols), struct2array(d_beta), 'mo')
title(strcat('Beta Band: R = ', num2str(R2(1,2)), ',', ' P = ', num2str(P2(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, ...
    '/home/javi/Documentos/measures_clustering_daniele/clustering/dist_vs_vols_beta.png')
close(f)

% Theta
[R3, P3] = corrcoef(struct2array(vols), struct2array(d_theta));

f = figure('visible','off');
scatter(struct2array(vols), struct2array(d_theta), 'rx')
title(strcat('Theta Band: R = ', num2str(R3(1,2)), ',', ' P = ', num2str(P3(1,2))))
xlabel('Volume (mm^3)')
ylabel('Average distance')
lsline
saveas(f, ...
    '/home/javi/Documentos/measures_clustering_daniele/clustering/dist_vs_vols_theta.png')
close(f)

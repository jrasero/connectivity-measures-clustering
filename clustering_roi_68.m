clc;clear;

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

    for iroi=1:n_rois
        for ihemi=1:2
            tic
            roi = rois_merged{iroi};
            hemi = hemispheres{ihemi};
            mkdir(strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/roi_68/',... 
                roi, "_", hemi), band)
            output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/clustering/roi_68/',...
                roi, "_", hemi, "/", band);
            R_rec = zeros(n_recs, n_measures*n_funcs, n_measures*n_funcs);

            for irec=1:n_recs
                load(strcat('similarity_matrices/roi_68', "/",...
                    roi,"_", hemi, "/", band, "/",  ...
                    'recording_', num2str(irec), '.mat'), 'R')
                  R_rec(irec, :, :) = R;
            end
            R_avg = squeeze(mean(R_rec, 1)); % average across recordings

             % 2-Plot this matrix
             f = figure('visible','off');
             imagesc(R_avg)
             xticks(1:n_measures*n_funcs)
             xticklabels(ticklabels)
             xtickangle(45)
             yticks(1:n_measures*n_funcs)
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
             toc
        end
    end
    toc
 end
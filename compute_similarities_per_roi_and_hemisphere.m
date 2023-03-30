clear; clc;
%% Define pipeline options

% Frequency bands
bands = {'alpha', 'beta', 'theta'}; 

% TS extraction methods
funcs = {'MaxAfter', 'PCABefore', 'MeanAfter', 'MeanBefore'}; 

% ROI labels
rois_merged = {'temporalpole', 'caudalmiddlefrontal', 'lateraloccipital', ...
    'inferiortemporal', 'supramarginal', ...
    'frontalpole', ...
    'rostralmiddlefrontal',...
    'middletemporal','fusiform',...
    'superiorparietal','pericalcarine', ...
    'superiortemporal',...
    'superiorfrontal',...
    'lingual',...
    'rostralanteriorcingulate',...
    'isthmuscingulate','entorhinal','parstriangularis',...
    'transversetemporal','bankssts',...
    'cuneus',...
    'posteriorcingulate','lateralorbitofrontal',...
    'precentral','parsopercularis',...
    'postcentral',...
    'caudalanteriorcingulate','inferiorparietal',...
    'parsorbitalis','insula','precuneus',...
    'paracentral',...
    'parahippocampal',...
    'medialorbitofrontal'
    };

% Connectivity measures
measures = {'PLV', 'ciPLI'};

% Hemispheres (Left/Righ)
hemispheres={'L', 'R'};

% Save these for future use
save('pipeline_options.mat', 'bands', 'funcs', 'rois_merged', 'measures', 'hemispheres')
%% Connectivity-based similarity matrix.
% Here we compute the similarity matrix between pipeline choices using
% the connectivity patterns per Region, concatenating their left and 
% right profiles.

n_bands = length(bands);
n_func = length(funcs);
n_measures = length(measures);
n_rois = length(rois_merged);
n_recs = 160;

% Iterate through (merged) regions
for iroi=1:n_rois
    for ihem=1:2
    tic
    roi = rois_merged{iroi};
    hemi = hemispheres{ihem};
    disp(['Doing ROI:' ' ' roi ', hemisphere: ' hemi])

    % Iterate through bands
    for iband=1:n_bands

      band = bands{iband};

      % Create output directory
      mkdir(strcat('/home/javi/Documentos/measures_clustering_daniele/similarity_matrices/roi_68/',... 
          [roi '_' hemi]), band)

      output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/similarity_matrices/roi_68/',...
           [roi '_' hemi], "/", band);

      % We compute a similarity matrix per recording. 
      % Later, prior to clustering we average across recordings.
      for irec=1:n_recs

        % Similarity matrix between analytical choices
        R = zeros(n_func*n_measures, n_func*n_measures);

        % Iterate through measure
        for ml=1:n_measures
            % Iterate through way of extracting the time series
            for fl=1:n_func
                % composite left side index
                ii = n_func*(ml-1)+fl;
                % Extract connectivity pattern by concatenating left
                % and hemisphere 
                folder = strcat(measures{ml}, "/", ...
                    measures{ml}, "_", roi, "_", hemi,"_", ...
                    band, "_", funcs{fl});
                files = dir(folder);
                load(strcat(files(2+irec).folder, "/", files(2+irec).name), 'TF')
                % load time series
                Cleft = transpose(TF);
                
                % Correlate this pattern with another anlytical option
                for mr=1:n_measures
                    for fr=1:n_func
                        % composite right side index
                        jj = n_func*(mr-1)+fr;
                        if (ii==jj)
                            R(ii,jj)=1;
                            continue
                        end
                        if(ii>jj)
                            continue
                        end


                        folder = strcat(measures{mr}, "/", ...
                            measures{mr}, "_", roi, "_",hemi,"_", ...
                            band, "_", funcs{fr});
                        files = dir(folder);
                        load(strcat(files(2+irec).folder, "/", files(2+irec).name), 'TF')
                        % Load time series
                        Cright = transpose(TF);

                        r = corrcoef(Cleft, Cright);
                        R(ii,jj) = r(1, 2);
                        R(jj,ii) = R(ii,jj);

                    end   
                end
            end
        end

        output_name = strcat(output_dir, '/', "recording_", ...
            num2str(irec), '.mat');
        save(output_name, 'R')
      end
      disp(['Band' ' ' band ' ' 'Done'])
     end
     toc
     disp(' ')
     
    end
end
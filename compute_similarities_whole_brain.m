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

% Data directory
data_dir = '/media/javi/ExtraDrive21/measures_clustering_daniele';

%% Connectivity-based similarity matrix.
% Here we compute the similarity matrix between pipeline choices using
% the connectivity patterns concatenated across seed regions.

n_bands = length(bands);
n_func = length(funcs);
n_measures = length(measures);
n_rois = length(rois_merged);
n_recs = 160;
  

% Iterate through bands
for iband=1:n_bands
  tic
  band = bands{iband};

  % Create output directory
  mkdir(strcat('/home/javi/Documentos/measures_clustering_daniele/similarity_matrices/',... 
      'whole_brain'), band)

  output_dir = strcat('/home/javi/Documentos/measures_clustering_daniele/similarity_matrices/',...
       'whole_brain', "/", band);

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
            Cleft = [];
            for iroi=1:n_rois
                roi = rois_merged{iroi};
                for ihem=1:2
                    folder = strcat(data_dir, "/", measures{ml}, "/", ...
                        measures{ml}, "_", roi, "_", hemispheres{ihem},"_", ...
                        band, "_", funcs{fl});
                    files = dir(folder);
                    load(strcat(files(2+irec).folder, "/", files(2+irec).name), 'TF')
                    % Concatenate data from both 
                    Cleft = [Cleft transpose(TF)];
                end
            end

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
                    Cright=[];
                    for iroi=1:n_rois
                        roi = rois_merged{iroi};
                        for ihem=1:2
                            folder = strcat(data_dir, "/", measures{mr}, "/", ...
                                measures{mr}, "_", roi, "_", hemispheres{ihem},"_", ...
                                band, "_", funcs{fr});
                            files = dir(folder);
                            load(strcat(files(2+irec).folder, "/", files(2+irec).name), 'TF')
                            % Concatenate data from both 
                            Cright = [Cright transpose(TF)];
                        end
                    end

                    r = corrcoef(Cleft, Cright);
                    R(ii,jj) = r(1,2);
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
  toc
end
%% Task-SRC
% The script is designed to analyze functional MRI (fMRI) data
% by correlating randomly permuted prediction models with neural activity patterns
% through representational similarity analysis (RSA).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permutations - Stelzer, Chen & Turner (2013) method
% Adapted for specific data analysis needs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Define the main directory, which is the common directory for all subject subfolders
home = 'E:/DataMRI/Scan';
cd(home);

% Define subject numbers
subj = [1:31, 33:36]; % There are 34 subjects
roisFILES = {'Fusiform_L', 'Fusiform_R', 'Putamen_L', 'Putamen_R'}; % Names of ROIs
nROIs = length(roisFILES); % Total number of ROIs

for r = 1:nROIs
    for S = 1:length(subj)
        % Enter each subject's specific model directory
        subj_dir = fullfile(home, ['subject' num2str(subj(S),'%.2d')...
            '/func/src/RSA/roiBased/funcROI/BasedRoiHundrVo/' roisFILES{r} '/']) ;
        cd(subj_dir);
        subj_cfg = [home, '/subject' num2str(subj(S),'%.2d')...
            '/func/src/RSA/roiBased/funcROI/BasedRoiHundrVo/' roisFILES{r} '/']; % Where should results be stored
        
        % Load the original configuration file
        cd(subj_cfg )
        load('res_cfg.mat'); % Assuming the original configuration file name is res_cfg.mat
        
        % Keep an unaltered cfg to copy parameters
        org_cfg = cfg;
        beta_loc = ['E:/DataMRI/Scan/subject' num2str(subj(S),'%.2d') '/func/src/FstLevSRC'];
        regressor_names = design_from_spm(beta_loc);
        loaded_data = load(fullfile(beta_loc,['residuals_' roisFILES{r} '.mat']));
        misc.residuals = loaded_data.misc.residuals;
        cfg.files.residuals.chunk = loaded_data.cfg.files.residuals.chunk;
        
        % Update cfg's results field
        cfg.results = rmfield(cfg.results, 'resultsname'); % Remove the results name field
        cfg.results.dir = fullfile(subj_cfg, 'perm/SRC_sum'); % Change the directory to the permutation directory
        if ~exist(cfg.results.dir)
            mkdir(cfg.results.dir)
        end
        ForSavePath = cfg.results.dir; % Path to save permutation results
        
        % Permutation design
        n_perms = 100; % Number of permutations
        files = org_cfg.files.name;
        files_perm = cell(n_perms, 1);
        for p = 1:n_perms
            sz = size(files, 1);
            rp = randperm(sz);
            files_perm{p, 1} = files(rp);
        end
        tic
        % Run all permutations in a loop
        for i_perm = 1:n_perms
            passed_data = []; % Avoid loading the same data multiple times
            disp(['Permutation ' num2str(i_perm) '/' num2str(n_perms)]);
            cfg.results.filestart = ['perm' sprintf('%04d', i_perm)];
            cfg.files.name = files_perm{i_perm, 1};
            % Perform decoding for this permutation:
            results = decoding(cfg,[],misc);
            % Rename the design files, starting with the current permutation number
            designfiles = dir(fullfile(cfg.results.dir, 'design.*'));
            for design_ind = 1:length(designfiles)
                movefile(fullfile(cfg.results.dir, designfiles(design_ind).name), ...
                    fullfile(cfg.results.dir, [cfg.results.filestart '_' designfiles(design_ind).name]));
            end
            Files{i_perm} = cfg.files.name;
        end
        toc
        % Clear variables and save permutation results
        clearvars -except home subj model S n_perms Files ForSavePath rois roisFILES r nROIs
        save(fullfile(ForSavePath, 'Permutations.mat'), 'Files');
    end
end

% Correlate the randomly permuted neural activity patterns with behavioral prediction models using RSA
% Specify the path of the behavioral prediction model
predictiveRDM_file = 'E:/DataMRI/CharRDMs/SRCsum.mat';
load(predictiveRDM_file); % Assuming this .mat file directly loads a variable named "SRCsum"

% Define the path for subjects and the file name pattern for permutation test results
subj_base_path = 'E:/DataMRI/Scan';
perm_file_pattern = 'perm%04d_other.mat';
n_permutations = 100; % Number of permutations
similarity_all = {};
for r = 1:nROIs
    
    % Process each subject in a loop
    for subj_num = 1:length(subj)
        % Use the dir function to get ROI file information
        tic
        subj_path = fullfile(subj_base_path, ['subject' sprintf('%.2d',subj(subj_num))...
            '/func/src/RSA/roiBased/funcROI/BasedRoiHundrVo/' roisFILES{r} '/']) ;
        % Process each permutation in a loop
        for perm_num = 1:n_permutations
            perm_file = fullfile(subj_path, 'perm/SRC_sum', sprintf(perm_file_pattern, perm_num));
            load(perm_file); % Load permutation test results
            
            % Calculate the similarity value for each voxel
            neural_activity = results.other.output{1, 1}; %{voxel_num, 1};
            
            % Create a lower triangular index matrix (excluding the diagonal)
            index_matrix = tril(ones(size(neural_activity)), -1);
            
            % Extract lower triangular values using the index matrix
            neural_activity = neural_activity(logical(index_matrix));
            rdm_pemutation_vector = SRCsum(logical(index_matrix));
            
            % Calculate Spearman correlation, and store the rho value
            rho = corr(neural_activity, rdm_pemutation_vector, 'Type', 'Spearman');
            
            % Perform Fisher-Z transformation
            z_val = 0.5 * log((1 + rho) / (1 - rho));
            similarity_values(perm_num) = z_val;
        end
        
        % Store the permutation test RSA distribution for 100 runs for the current subject and current ROI
        similarity_all{subj(subj_num),r+1} = similarity_values;
        similarity_all{subj(subj_num),1} = sprintf('subj%02d', subj(subj_num));
        disp([sprintf('subject%02d', subj(subj_num)),',done!' ])
        toc
    end
end

% Save the permutation test distribution results for all subjects and all ROIs
% Add column names
columnNames = {'Fusiform_L', 'Fusiform_R', 'Putamen_L', 'Putamen_R'};
% Create table
resultsTable = array2table(similarity_all, 'VariableNames', columnNames);
% Save the table
save('E:/DataMRI/RSA/MeanROIbased/SRC/BasedRoiHundrVo_SRCsum_Perm_Results.mat', 'resultsTable');



% This script runs the 2nd level (across participants) of the permutation procedure
% from each subject's each voxel's 100 random similarity values, one is randomly chosen and their mean is taken as the group level random similarity value result.
% this process is repeated 10,000 times, obtaining 10,000 group level random similarity value results. This is used to determine the p-value of the true value.
% then, FWE correction is performed. Results are saved in two files: one contains 10,000 group level random similarity values for all voxels (.mat format);
% the other contains each voxel's mask-index, real observed values, p-values (before FWE correction), and p-values (after FWE correction).
%% Fusiform_L SRC_left
clear; clc;

% Parameter settings
home = 'E:/DataMRI/Scan'; % Main directory
subj = [1:31, 33:36]; % Subject IDs
nSubj = length(subj); % Total number of subjects
% Load the predictive model RDM
predictiveRDM_file = 'E:/DataMRI/CharRDMs/SRC_left.mat';
load(predictiveRDM_file, 'SRC_left'); % This .mat file directly loads a variable named "SRC_left"
% Load the actual similarity coefficients and corresponding voxel index
load('E:/DataMRI/RSA/MeanROI_searchlight/SRC/Fusiform_L/BasedRoiHundrVo_Mean_SRC_left_RSA_zscores.mat');
observedVoxel = data;
nPermsFirstLevel = 100; % Number of first-level permutations
nPermsSecondLevel = 10000; % Number of second-level permutations

% Initialize storage for each second-level permutation group mean
groupPermResults = zeros(length(observedVoxel.mask_index), nPermsSecondLevel);
collect_vox = observedVoxel.mask_index;

% Second-level permutation loop
for Voxel = 1:length(observedVoxel.mask_index)
    currVox = observedVoxel.mask_index(Voxel);
    perm_zs = zeros(nSubj, nPermsSecondLevel);

    for perm2 = 1:nPermsSecondLevel
        tic
        % Loop over each subject
        for s = 1:nSubj
            subjID = subj(s);
            Subj_Perms = dir(fullfile(home, sprintf('subject%02d', subjID),...
                'func', 'src', 'RSA', 'ROI_searchlight', 'BasedRoiHundrVo',...
                'Fusiform_L', 'perm', 'SRC_sum', '*_other.mat'));
            
            % Randomly select one of the first-level permutation results for this subject
            randPermIndex = randi(nPermsFirstLevel);
            randPermRes = load(fullfile(Subj_Perms(randPermIndex).folder, ...
                Subj_Perms(randPermIndex).name));
            
            index = find(randPermRes.results.mask_index == currVox); % Find the current voxel's cell for this subject
            
            if ~isempty(index)
                permRDM = randPermRes.results.other.output{index, 1};
                % Convert the permutation activity matrix to a vector
                perm_vector = permRDM(tril(true(size(permRDM)), -1));
                
                % Calculate Spearman correlation with the predictive model
                perm_rho = corr(perm_vector, SRC_left(tril(true(size(SRC_left)), -1)), 'Type', 'Spearman');
                
                % Convert to Fisher Z-score
                perm_zs(s, perm2) = 0.5 * log((1 + perm_rho) / (1 - perm_rho));
            end
        end
        
        % Calculate and store the group mean for the current Voxel permutation
        groupPermResults(Voxel, perm2) = mean(perm_zs(:, perm2));
        
        disp( ['perm_' num2str(perm2) ', done!'])
    end
        disp(['NO.' num2str(Voxel) ', Vox_' num2str(currVox) ', done!'])
        toc

end

% All subjects all voxels' 10,000 permutation results
GroupPermVoxRes.PermVoxRes = groupPermResults;
GroupPermVoxRes.collect_vox = collect_vox; % Store voxel-index

% Save group level random similarity value results
save('E:/DataMRI/RSA/MeanROI_searchlight/SRC/Fusiform_L_SRC_left_GroupPermVoxRes.mat', 'GroupPermVoxRes');

% FWE correction
% Save files containing each voxel's mask-index, real observed values, p-values (before FWE correction), and p-values (after FWE correction)
% Observed Z-scores are stored in the observedVoxel structure's z_scores field
observedZ = observedVoxel.spearman_z_values;

% Initialize p-values array
pValues = zeros(size(observedZ));

% Calculate p-value for each voxel
for Voxel = 1:length(observedZ)
    pValues(Voxel) = round(mean(GroupPermVoxRes.PermVoxRes(Voxel, :) >= observedZ(Voxel)), 6);
end

% FWE correction
nPerms = size(GroupPermVoxRes.PermVoxRes, 2); % Total number of permutations
maxStatsPerPerm = max(GroupPermVoxRes.PermVoxRes, [], 2); % Find the maximum statistical value for each permutation

% For each Voxel, calculate the proportion of maximum statistical values greater than the observed mean, i.e., the p-value after FWE correction
fweCorrectedPValues = zeros(size(observedZ));
for Voxel = 1:length(observedZ)
    fweCorrectedPValues(Voxel) = round(mean(maxStatsPerPerm > observedZ(Voxel)),32);
end

% Save relevant data
for Voxel = 1:length(observedZ)
    % Check if observedZ is NaN
    if ~isnan(observedZ(Voxel))
        % Build table to save data
        T(Voxel,:) = table(collect_vox(Voxel), observedZ(Voxel), pValues(Voxel), fweCorrectedPValues(Voxel),...
            'VariableNames', {'Mask_Index', 'Observed_Z', 'P_Value', 'fweCorrectedPValues'});
    end
end
    % Save data to CSV file
    writetable(T, fullfile('E:/DataMRI/RSA/MeanROI_searchlight/SRC/Fusiform_L_SRC_left_GroupPermVoxRes_fweCorr.csv'));


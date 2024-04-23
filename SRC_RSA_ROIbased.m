%% Task-SRC RSA
% This script is a template that can be used for a representational
% similarity analysis on brain image data. It is for people who have betas
% available from an SPM.mat and want to automatically extract the relevant
% images used for the similarity analysis, as well as corresponding labels
% and decoding chunk numbers (e.g. run numbers). If you don't have this
% available, then inspect the differences between decoding_template and
% decoding_template_nobetas and adapt this template to use it without
% betas.

clear all;
clc;

% set group parameters
subdir = [1:31 33:36]; % subject 32 is excluded due to excessive head movement
filebase = 'E:/DataMRI/Scan/subject'; % common folder where all subject subfolders are located (e.g., 'D:/data/sub')
corr_type = 'Pearson';% distance measure (usually Pearson's r -- to see other distance options do 'edit pattern_similarity' in the command window and check lines 14:23)
roisFILES = {'Fusiform_L', 'Fusiform_R', 'Putamen_L', 'Putamen_R'}; % Names of ROIs
nROIs = length(roisFILES); % Total number of ROIs

for r = 1:nROIs
    
    for i = 1:length(subdir)
        
        % Set defaults
        cfg = decoding_defaults;
        cfg.plot_design=0;
        cfg.results.overwrite = 1;
        cfg.analysis = 'ROI';
        fileroot = ([filebase num2str(subdir(1,i),'%.2d')  '/']);
        
        % start looping over rois
        cfg.results.dir = [filebase num2str(subdir(1,i),'%.2d')...
            '/func/src/RSA/roiBased/funcROI/BasedRoiHundrVo/' roisFILES{r} '/'] ;  % where should results be stored
        
        beta_loc = [fileroot 'func/src/FstLevSRC']; % beta images should be here
        cfg.files.mask = ([filebase num2str(subdir(1,i),'%.2d')...
            '/func/src/BasedRoiHundrVoFst_GLM/' roisFILES{r} '/hundrVoxel.nii']);  % mask = roi in this case, rois are resliced to participants native space.
        regressor_names = design_from_spm(beta_loc);
        
        % If you didn't specifiy a label order before, set the label names to
        % the regressor names which you want to use here, e.g. 'button left' and 'button right'
        % don't remember the names? -> run display_regressor_names(beta_loc)
        labelnames = {'SRCtarg001', 'SRCtarg002', 'SRCtarg003', 'SRCtarg004', 'SRCtarg005', 'SRCtarg006', 'SRCtarg007', 'SRCtarg008', 'SRCtarg009', 'SRCtarg010',...
            'SRCtarg011', 'SRCtarg012', 'SRCtarg013', 'SRCtarg014', 'SRCtarg015', 'SRCtarg016', 'SRCtarg017', 'SRCtarg018', 'SRCtarg019', 'SRCtarg020', 'SRCtarg021', ...
            'SRCtarg022', 'SRCtarg023', 'SRCtarg024', 'SRCtarg025', 'SRCtarg026', 'SRCtarg027', 'SRCtarg028', 'SRCtarg029', 'SRCtarg030', 'SRCtarg031','SRCtarg032',...
            'SRCtarg033', 'SRCtarg034', 'SRCtarg035', 'SRCtarg036'};
        
        % make con-images in beta_loc. Each trial repeated six times, and we averaged them to get the mean beta map for each trial.
        for i_label = 1:length(labelnames)
            % get corresponding beta from regressor_names.mat
            beta_ind = find(contains(regressor_names(:),labelnames{i_label}));
            if length(beta_ind) < 6
                error('We need 6 matching betas for condition %s but only %i found. Please check the naming.',labelnames{i_label},length(beta_ind))
            end
            hdr1 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(1))]);
            hdr2 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(2))]);
            hdr3 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(3))]);
            hdr4 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(4))]);
            hdr5 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(5))]);
            hdr6 = spm_vol([beta_loc '/' sprintf('beta_%04i.nii',beta_ind(6))]);
            
            vol1 = spm_read_vols(hdr1);
            vol2 = spm_read_vols(hdr2);
            vol3 = spm_read_vols(hdr3);
            vol4 = spm_read_vols(hdr4);
            vol5 = spm_read_vols(hdr5);
            vol6 = spm_read_vols(hdr6);
            
            vol7 = 1/6 * (vol1+vol2+vol3+vol4+vol5+vol6);
            hdr7 = hdr6;
            hdr7.fname = fullfile(beta_loc,sprintf('con_%04i.nii',i_label));
            spm_write_vol(hdr7,vol7)
        end
        
        % now set up analysis (this replaces decoding_describe_data)
        cfg.files.name = {};
        for i_label = 1:36
            cfg.files.name{i_label,1} = fullfile(beta_loc,sprintf('con_%04i.nii',i_label));
        end
        cfg.files.label = (2*mod(1:36,2)-1)'; % just 1 and -1 since it is arbitrary
        cfg.files.chunk = ones(36,1);
        
        % set everything to similarity analysis
        cfg.decoding.software = 'similarity';
        cfg.decoding.method = 'classification';
        cfg.decoding.train.classification.model_parameters = corr_type;
        cfg.results.output = 'other';
        
        % set normalization
        % These parameters carry out the multivariate noise normalization using the
        % residuals.
        % The crossnobis distance is identical to the cross-validated Euclidean
        % distance after prewhitening (multivariate noise normalization). It has
        % been shown that a good estimate for the multivariate noise is provided
        % by the residuals of the first-level model, in addition with Ledoit-Wolf
        % regularization. Here we calculate those residuals. If you have them
        % available already, you can load them into misc.residuals using only the
        % voxels from cfg.files.mask
        if ~exist(fullfile(beta_loc,['residuals_' roisFILES{r} '.mat']),'file')
            cfg.scale.method = 'cov'; % we scale by noise covariance
            cfg.scale.estimation = 'separate'; % we scale all data for each run separately while iterating across searchlight spheres
            cfg.scale.shrinkage = 'lw2'; % Ledoit-Wolf shrinkage retaining variances
            [misc.residuals,cfg.files.residuals.chunk] = residuals_from_spm(fullfile(beta_loc,'SPM.mat'),cfg.files.mask); % this only needs to be run once and can be saved and loaded
            save((fullfile(beta_loc,['residuals_' roisFILES{r} '.mat'])),'misc', 'cfg')
        else
            loaded_data = load(fullfile(beta_loc,['residuals_' roisFILES{r} '.mat']));
            misc.residuals = loaded_data.misc.residuals;
            cfg.files.residuals.chunk = loaded_data.cfg.files.residuals.chunk;
        end
        cfg.verbose = 0;
        
        %% Nothing needs to be changed below for a standard similarity analysis using all data
        cfg.design = make_design_similarity(cfg);
        cfg.design.unbalanced_data = 'ok';
        
        % Run decoding
        results = decoding(cfg,[],misc);
        
        % Save the current ROI results to the all_results array
        all_ROI_results(r).ROI = roisFILES{r};
        all_ROI_results(r).data = results;
    end
    % Save all_results to a .mat file
    save([filebase num2str(subdir(1,i),'%.2d') '/func/src/RSA/roiBased/funcROI/' 'BasedRoiHundrVo_all_ROI_results.mat'], 'all_ROI_results');
end


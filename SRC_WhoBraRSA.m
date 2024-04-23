%% Task-SRC RSA searchlight based on whole brain
% This script is a template that can be used for a representational
% similarity analysis on brain image data. It is for people who have betas
% available from an SPM.mat and want to automatically extract the relevant
% images used for the similarity analysis, as well as corresponding labels
% and decoding chunk numbers (e.g. run numbers). If you don't have this
% available, then inspect the differences between decoding_template and
% decoding_template_nobetas and adapt this template to use it without
% betas.

clear all
clc
subdir = [1:31 33:36]; % subject 32 is excluded due to excessive head movement
filebase =  'E:/DataMRI/Scan/subject'; % common folder where all your subjects subfolders are (e.g. 'D:/data/sub')
corr_type = 'pearson'; % distance measure (usually Pearson's r -- to see other distance options do 'edit pattern_similarity' in the command window and check lines 14:23)

for i = 1 : length(subdir)
    tic
    % Set defaults
    cfg = decoding_defaults;
    cfg.plot_design=0;
    cfg.results.overwrite = 1;
    cfg.analysis = 'searchlight';
    cfg.searchlight.unit = 'voxels';
    cfg.searchlight.radius = 6;
    fileroot = [filebase sprintf('%02d/func/src', subdir(i))];
    cfg.results.dir = ([fileroot '/RSA/WholeBrain' ]); % where should results be stored
    beta_loc = [fileroot '/FstLevSRC']; % beta images should be here
    cfg.files.mask = strcat(filebase,sprintf('%02d/func/src/FstLevSRC', subdir(i)), "/mask.nii"); % 'E:/DataMRI/ROI/rEPI.nii';  this is the whole brain (output of SPM GLM)
    load(fullfile(beta_loc,'regressor_names.mat'))
    
    labelnames =  {'SRCtarg001', 'SRCtarg002', 'SRCtarg003', 'SRCtarg004', 'SRCtarg005', 'SRCtarg006', 'SRCtarg007', 'SRCtarg008', 'SRCtarg009', 'SRCtarg010',...
        'SRCtarg011', 'SRCtarg012', 'SRCtarg013', 'SRCtarg014', 'SRCtarg015', 'SRCtarg016', 'SRCtarg017', 'SRCtarg018', 'SRCtarg019', 'SRCtarg020', 'SRCtarg021', ...
        'SRCtarg022', 'SRCtarg023', 'SRCtarg024', 'SRCtarg025', 'SRCtarg026', 'SRCtarg027', 'SRCtarg028', 'SRCtarg029', 'SRCtarg030', 'SRCtarg031','SRCtarg032',...
        'SRCtarg033', 'SRCtarg034', 'SRCtarg035', 'SRCtarg036'};
    
    % make con-images in beta_loc. Each trial repeated twice, and we averaged them to get the mean beta map for each trial.
    for i_label = 1:length(labelnames)
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
    % Decide whether you want to see the searchlight/ROI/... during decoding
    cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...
    cfg.verbose = 0; % change to 1 or 2 if you want to get more feedback while the script is running
    
    %% Nothing needs to be changed below for a standard similarity analysis using all data
    regressor_names = design_from_spm(beta_loc);
    cfg.design = make_design_similarity(cfg);
    cfg.design.unbalanced_data = 'ok';
    
    % Run decoding
    results = decoding(cfg);
    toc
    
    disp(strcat('subj0 ',string(i),' finished'))
end


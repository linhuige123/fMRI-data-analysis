%% Task-SRC Preprocessing
% This script is designed for preprocessing Magnetic Resonance Imaging (MRI) data.
% Data Reading and Path Configuration
% Iterating Over Subjects: The script is designed for 36 subjects, iterating through each subject.
% NIfTI File Reading: Retrieves paths for NIfTI files (`.nii`) in each run session directory, which contain raw scan data.
% JSON Configuration Reading: Configuration files for each run (like scan timing parameters) are stored in JSON files, which the script parses to obtain scanning parameters.

% SPM Batch Configuration
% Set Up Batch Tasks: Configures image processing batch tasks using the MATLAB interface of SPM (Statistical Parametric Mapping) software.
% Slice Timing Correction: Applies timing corrections using parameters from the JSON file.
% Motion Correction: Corrects images for motion to ensure image quality.
% Image Coregistration: Coregisters functional images to anatomical images.
% Tissue Segmentation and Normalization: Segments anatomical images into different tissue types and normalizes images to a standard brain template.

%-----------------------------------------------------------------------
% Job saved on 06-Jan-2023 17:29:27 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

clear;
cwd = 'E:/DataMRI/Scan';

for subject = 1:36 %subjects
    
    if subject<10
        CurrSub = [cwd '/' sprintf('subject%02d/',subject) ]; %  path of the current subject
        subject = num2str(subject, '%02d');
    else
        CurrSub = [cwd '/' sprintf('subject%01d/',subject) ]; %  path of the current subject
        subject = num2str(subject, '%01d');
    end

Run = dir([CurrSub '/func/opc/Proprecessing/']) ;
Run = Run(3:end);
Run1 = strcat(Run(1).folder, '/run1'); 
Run2 = strcat(Run(1).folder, '/run2'); 
Run3 = strcat(Run(1).folder, '/run3'); 
Run4 = strcat(Run(1).folder, '/run4'); 
Run5 = strcat(Run(1).folder, '/run5'); 
Run6 = strcat(Run(1).folder, '/run6'); 

fileDir1 = dir([Run1 '/'  'run*.nii']); file1 = [fileDir1.folder '/' fileDir1.name ];
fileDir2 = dir([Run2 '/' 'run*.nii']); file2 = [fileDir2.folder '/' fileDir2.name ];
fileDir3 = dir([Run3 '/' 'run*.nii']); file3 = [fileDir3.folder '/' fileDir3.name ];
fileDir4 = dir([Run4 '/' 'run*.nii']); file4 = [fileDir4.folder '/' fileDir4.name ];
fileDir5 = dir([Run5 '/' 'run*.nii']); file5 = [fileDir5.folder '/' fileDir5.name ];
fileDir6 = dir([Run6 '/' 'run*.nii']); file6 = [fileDir6.folder '/' fileDir6.name ];

json1=dir([fileDir1.folder '/*.JSON']);
json=jsondecode(fileread(fullfile(json1.folder,json1.name)));
[~,pos]=min(abs(json.SliceTiming-0.5));

T1 = dir([CurrSub '/anat/T1/']) ; T1Path = strcat(T1(1).folder);
T1fileDir = dir([T1Path '/T1' '*.nii']); T1file = [T1fileDir.folder '/' T1fileDir.name ];

% matlabbatch
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'opcRun1Run6';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                     {file1}
                                                                     {file2}
                                                                     {file3}
                                                                     {file4}
                                                                     {file5}
                                                                     {file6}
                                                                     }';
matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Named File Selector: opcRun1Run6(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{2}.spm.temporal.st.scans{2}(1) = cfg_dep('Named File Selector: opcRun1Run6(2) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
matlabbatch{2}.spm.temporal.st.scans{3}(1) = cfg_dep('Named File Selector: opcRun1Run6(3) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{3}));
matlabbatch{2}.spm.temporal.st.scans{4}(1) = cfg_dep('Named File Selector: opcRun1Run6(4) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{4}));
matlabbatch{2}.spm.temporal.st.scans{5}(1) = cfg_dep('Named File Selector: opcRun1Run6(5) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{5}));
matlabbatch{2}.spm.temporal.st.scans{6}(1) = cfg_dep('Named File Selector: opcRun1Run6(6) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{6}));
matlabbatch{2}.spm.temporal.st.nslices = length(json.SliceTiming);
matlabbatch{2}.spm.temporal.st.tr = json.RepetitionTime;
matlabbatch{2}.spm.temporal.st.ta = 0;
matlabbatch{2}.spm.temporal.st.so = json.SliceTiming*1000;
matlabbatch{2}.spm.temporal.st.refslice = json.SliceTiming(pos(1))*1000;
matlabbatch{2}.spm.temporal.st.prefix = 'a';
matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.data{3}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 3)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{3}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.data{4}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 4)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{4}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.data{5}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 5)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{5}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.data{6}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{6}, '.','files'));
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
matlabbatch{4}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{4}.spm.spatial.coreg.estimate.source = {T1file};
matlabbatch{4}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{5}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{5}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{5}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{5}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{5}.spm.spatial.preproc.tissue(1).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,1'};
matlabbatch{5}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{5}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{5}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(2).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,2'};
matlabbatch{5}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{5}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{5}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(3).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,3'};
matlabbatch{5}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{5}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{5}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(4).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,4'};
matlabbatch{5}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{5}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{5}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(5).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,5'};
matlabbatch{5}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{5}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{5}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(6).tpm = {'C:/Program Files/MATLAB/R2020a/toolbox/spm12/tpm/TPM.nii,6'};
matlabbatch{5}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{5}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{5}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{5}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{5}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{5}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{5}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{5}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{5}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{5}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{5}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{5}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{6}.spm.spatial.normalise.write.subj(1).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(1).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.subj(2).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(2).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.subj(3).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(3).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.subj(4).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(4).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.subj(5).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(5).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 5)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.subj(6).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{6}.spm.spatial.normalise.write.subj(6).resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 6)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','cfiles'));
matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w';

% data for RSA need not to be smoothed
% matlabbatch{7}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.data(2) = cfg_dep('Normalise: Write: Normalised Images (Subj 2)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.data(3) = cfg_dep('Normalise: Write: Normalised Images (Subj 3)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.data(4) = cfg_dep('Normalise: Write: Normalised Images (Subj 4)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.data(5) = cfg_dep('Normalise: Write: Normalised Images (Subj 5)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.data(6) = cfg_dep('Normalise: Write: Normalised Images (Subj 6)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{7}.spm.spatial.smooth.fwhm = [6 6 6];
% matlabbatch{7}.spm.spatial.smooth.dtype = 0;
% matlabbatch{7}.spm.spatial.smooth.im = 0;
% matlabbatch{7}.spm.spatial.smooth.prefix = 's';

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);
end

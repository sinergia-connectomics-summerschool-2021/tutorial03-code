clc
close all
clear

% make sure you are in the correct folder
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))

addpath(genpath('.'))

try
    ft_defaults % check if fieldtrip is in Matlab path
catch
    addpath(fullfile('/','home','sinergiasummerschool','Softwares','fieldtrip'))
    ft_defaults
end

%% initialisation
%subj init
BIDS_folder=fullfile('/', 'home','sinergiasummerschool','Data','ds003505');
MRI_preproc_folder=fullfile(BIDS_folder,'derivatives','cmp');
EEG_preproc_folder=fullfile(BIDS_folder,'derivatives','eeg_preprocessing');
task='task-faces';
id_sub=1;
sub_id = sprintf('sub-%02d',id_sub);

%mri file
mri_filename = fullfile(MRI_preproc_folder, sub_id,'anat', [sub_id,'_desc-head_T1w.nii.gz']); 
%gray matter map file
GM_filename=fullfile(EEG_preproc_folder,sub_id,'anat',[sub_id,'_label-GM_dseg.nii.gz']); 
%atlas file
atlas_filename=fullfile(MRI_preproc_folder,sub_id,'anat',[sub_id,'_atlas-L2018_res-scale1_dseg.nii.gz']);  
%electrode file (coregistered to mri) 
elec_file=fullfile(BIDS_folder,sub_id, 'eeg', [sub_id '_electrodes.tsv']);

%preproc init
%headmodel
cfg_headmdl.tissues={'brain','skull','scalp'}; %band stop freq
cfg_headmdl.numvertices=[3000,2000,1000]; %vertices for tissues
cfg_headmdl.method='bemcp';
%sourcemodel
cfg_sourcemdl.ROI2remove=[35:39,76:80, 83]; %thalamus, caudate, putamen, pallidum, accumbens area, brainstem
cfg_sourcemdl.method_reslice='cubic';
cfg_sourcemdl.resolution_reslice=1; %1 x 1 x 1 mm
cfg_sourcemdl.grid_resolution=6; %source grid resolution (6mm)
cfg_sourcemdl.unit='mm'; %unit of source grid

%% CREATE HEADMODEL: segmentation, mesh, volume conduction model
%read mri and smooth it to the segmentation
[mri_fname, temp_folder] = decompress_mri(mri_filename);
mri = ft_read_mri(mri_fname);
mri.coordsys = 'ras'; % we suppose that the coord sys is ras
delete(mri_fname);    rmdir(temp_folder);

%-segment mri
cfg           = [];
cfg.output = {'brain', 'skull'};
segmentedmri  = ft_volumesegment(cfg, mri);

% first little trick applied here because the mri has been defaced, it makes sure the scalp has no holes
cfg           = [];
cfg.output = {'scalp'};
segmentedmri2  = ft_volumesegment(cfg, mri);

segmentedmri.scalp = segmentedmri2.scalp;
clear segmentedmri2

segmentedmri.scalp = segmentedmri.scalp | ...
    imdilate(segmentedmri.skull, strel('sphere',6)); % 2nd little trick applied here to recreate a bit of the forehead

%create surfaces (mesh at tissues interfaces)
cfg=[];
cfg.tissue=cfg_headmdl.tissues;
cfg.numvertices=cfg_headmdl.numvertices;
mesh=ft_prepare_mesh(cfg,segmentedmri);

%create a volume conduction model
cfg=[];
cfg.method=cfg_headmdl.method;
headmodel  = ft_prepare_headmodel(cfg, mesh);

%% CREATE SOURCEMODEL: remove basal ganglia from GM, coregister to tmeplate and reslice, impose sources inside GM
% Read gray matter MRI file
[gray_fname, temp_folder] = decompress_mri(GM_filename);
gray = ft_read_mri(gray_fname);
gray.coordsys = 'ras'; % we suppose that the coord sys is ras
delete(gray_fname);     rmdir(temp_folder);

% Read ATLAS file
[atlas_fname, temp_folder] = decompress_mri(atlas_filename);
atlas = ft_read_mri(atlas_fname);
atlas.coordsys = 'ras'; % we suppose that the coord sys is ras
delete(atlas_fname); rmdir(temp_folder);

% mask GM map with ROIs2remove from GM
mask=ismember(atlas.anatomy,cfg_sourcemdl.ROI2remove);
gray.anatomy(mask)=0;

[~, ftpath] = ft_version;

% Read template MRI file
tmp_mri_filename = fullfile(ftpath, 'external', 'spm12','toolbox', 'OldNorm', 'T1.nii');
tmp_mri = ft_read_mri(tmp_mri_filename);
tmp_mri.coordsys = 'ras';

% Co-register the mri to the template
cfg = [];
cfg.method = 'spm';
cfg.spm.regtype = 'rigid';
cfg.viewresult = 'no';
[realign] = ft_volumerealign(cfg, mri, tmp_mri);

% set the same transformation to the gray matter prob map
realign_gray = gray;
realign_gray.transform = realign.transform;

% apply 1 mm reslicing
cfg=struct('method',cfg_sourcemdl.method_reslice,'resolution',cfg_sourcemdl.resolution_reslice,'dim',250*[1 1 1]);
realign_reslice_GM = ft_volumereslice(cfg, realign_gray);
clear realign_gray

% create source model restricted to gray matter
% sourcemodel = prepare_sourcemodel(realign_reslice_GM, inputs);
cfg=struct('resolution',cfg_sourcemdl.grid_resolution,'unit',cfg_sourcemdl.unit,'tight','yes','mri',realign_reslice_GM);
sourcemodel = ft_prepare_sourcemodel(cfg);

%count sources in the brain
num_sources_in_brain = sum(sourcemodel.inside);

% transform the coordinate of the source model back to the original mri
% space (they are currently in the space of the template)
T_mat = realign.transformorig/realign.transform;
[sourcemodel.pos]= ft_warp_apply(T_mat, sourcemodel.pos, 'homogeneous');

%% DISPLAY RESULTS
elec = ft_read_sens(elec_file); %yet to be coregistered
fig = verify_head_and_source_models(headmodel, sourcemodel, elec, 1);

%% save in BIDS

derivatives_path = fullfile(BIDS_folder, 'derivatives', 'mri_preprocessing', sub_id, 'anat');
if ~exist(derivatives_path, 'dir')
   mkdir(derivatives_path) 
else
    % Tell Datalad to allow files to be modified
    mri_preprocessing_derivatives_dir = fullfile(BIDS_folder, 'derivatives', 'mri_preprocessing');
    [status,cmdout] = system('datalad unlock -d '+convertCharsToStrings(BIDS_folder)+' '+ convertCharsToStrings(mri_preprocessing_derivatives_dir));
    sprintf(cmdout)
end

%save headmodel
headmodel_fname=[sub_id,'_space-individual_desc-reslice_headmodel'];
save(fullfile(derivatives_path,headmodel_fname),'headmodel');

%save sourcemodel
sourcemodel_fname=[sub_id,'_space-individual_desc-reslice_sourcemodel'];
save(fullfile(derivatives_path,sourcemodel_fname),'sourcemodel');

% store image of the results
photo_fname = [sub_id,'_space-individual_desc-reslice_photo.jpg'];
saveas(fig, fullfile(derivatives_path,photo_fname));



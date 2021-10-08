%This script is used to :
%   - load eeg data, headmodel and sourcemodel 
%   - get covariance matrix on baseline eeg
%   - compute leadfield (forward solution)
%   - estimate the inverse filters with MNE (inverse solution)

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
BIDS_folder = fullfile('/', 'home','sinergiasummerschool','Data','ds003505');
derivatives_folder=fullfile('/', 'home','sinergiasummerschool','Data','ds003505', 'derivatives');
task='task-faces';
id_sub = 1;
sub_id = sprintf('sub-%02d',id_sub);
eeg_fname = fullfile(derivatives_folder, 'eeg_preprocessing', sub_id, 'eeg','data'); %eeg file
headmodel_fname=fullfile(derivatives_folder, 'mri_preprocessing',sub_id,'anat', [sub_id,'_space-individual_desc-reslice_headmodel']); %headmodel file
sourcemodel_fname=fullfile(derivatives_folder, 'mri_preprocessing',sub_id,'anat',[sub_id,'_space-individual_desc-reslice_sourcemodel']); %sourcemodel file
load('elec.mat')

derivatives_path = fullfile(derivatives_folder, 'source_modelling', sub_id, 'eeg');
if ~exist(derivatives_path, 'dir')
   mkdir(derivatives_path) 
end


%preproc init 
cfg_inverse.lambda =.05; 
cfg_inverse.method = 'eloreta';

%% load stuff
%headmodel
load(headmodel_fname)
%sourcemodel
load(sourcemodel_fname)
%preprocessed eeg
load(eeg_fname)

%% Tell Datalad to allow files to be modified
if exist(derivatives_path, 'dir')
    [status,cmdout] = system('datalad unlock -d '+convertCharsToStrings(BIDS_folder)+' '+ convertCharsToStrings(derivatives_path));
    sprintf(cmdout)
end

%% estimate the leadfield
cfg         = [];
cfg.elec    = elec_proj;   % sensor information
cfg.sourcemodel    = sourcemodel;   % source points
cfg.headmodel = headmodel;   % volume conduction model
leadfield   = ft_prepare_leadfield(cfg);

%save leadfield
leadfield_fname=[sub_id,'_space-individual_desc-cleaned_leadfield'];
save(fullfile(derivatives_path,leadfield_fname),'leadfield');
%% estimate the inverse operator
dat = []; % we only want the filter for now
headmodel = []; % the headmodel is only used when isinside not present in leadfield and if sourcemodel does not have leadfield
elec = []; % the elec structure is only need if leadfield not present which is not hte case here

inverse_operator = ft_inverse_eloreta(leadfield, elec, headmodel, dat, ...
    eye(length(data.label)), 'keepfilter', 'yes',...
    'lambda', cfg_inverse.lambda);

%save inverse filters
filter_fname=[sub_id ,'_space-individual_desc-cleaned_inverseoperator'];
save(fullfile(derivatives_path,filter_fname),'inverse_operator');

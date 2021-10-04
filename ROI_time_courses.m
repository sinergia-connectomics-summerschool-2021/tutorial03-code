%This script is used to project the EEG data in source space and create the ROIs traces:
%   - import the cleaned EEG data 
%   - import the electrode position
%   - import the inverse filter
%   - import the atlas in the mri space
%   - import table with ROI names
%   - project the data in source space
%   - create the ROI time-courses
%   - store the ROI time-courses

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
% subj init
BIDS_folder=fullfile('/', 'home','sinergiasummerschool','Data','ds003505');

MRI_preproc_folder=fullfile(BIDS_folder,'derivatives', 'cmp');
task='task-faces';
id_sub=1;
sub_name = sprintf('sub-%02d',id_sub);
clean_eeg_list = dir(fullfile(BIDS_folder,'derivatives','eeg_preprocessing',...
    sub_name, '**', 'data.mat'));
clean_eeg_fname = fullfile(clean_eeg_list.folder,clean_eeg_list.name);

atlas_list = dir(fullfile(BIDS_folder,'derivatives','cmp',...
    sub_name, '**', [sub_name '_atlas-L2018*scale1_dseg.nii.gz']));
atlas_fname = fullfile(atlas_list.folder,atlas_list.name);

invop_list = dir(fullfile(BIDS_folder,'derivatives','source_modelling',...
    sub_name, '**', '*inverseoperator.mat'));
invfilters_fname = fullfile(invop_list.folder,invop_list.name);

clear clean_eeg_list atlas_list filters_list

ROI2remove=[35:39,76:80, 83];

derivatives_path = fullfile(BIDS_folder,'derivatives', 'source_modelling', sub_name, 'eeg');

%% import all data

%read atlas
[unzip_atlas_fname, temp_folder] = decompress_mri(atlas_fname);
atlas = ft_read_mri(unzip_atlas_fname);
delete(unzip_atlas_fname); 
rmdir(temp_folder);
clear unzip_atlas_fname temp_folder

%load inverse filters
load(invfilters_fname, 'inverse_operator');

%load clean EEG (it contains the electrode positions)
load(clean_eeg_fname, 'data');
%% Read table with ROIs and assign label to point source
Tbl_idx2label = readtable('Lausanne_atlas_label.tsv','FileType',...
            'text', 'Delimiter', '\t');

mask=ismember(atlas.anatomy,ROI2remove);
atlas.anatomy(mask) = 0;

[src_label, ROIs_elec]= ...
    assigned_label2source_points(inverse_operator, atlas, Tbl_idx2label);

%% 
ROIs_lbl = unique(src_label(~strcmp(src_label,'Outside')), 'stable');
[bln, idx] = ismember(Tbl_idx2label.abbreviation,ROIs_lbl);
ROIs_lbl = ROIs_lbl(nonzeros(idx));

n_ROIs = length(ROIs_lbl);
filters = inverse_operator.filter;
Fs = data.fsample;

[n_channels, n_samples] = size(data.trial{1}); % only works if trials have the same length
n_trials = length(data.trial);
time_courses = cat(2, data.trial{:});

ROI_traces = zeros(n_ROIs, n_samples*n_trials);

% parfor uses the parallel toolbox to parallelise the loop, change to for if this toolbox is not available
parfor ROI_id = 1:n_ROIs
    curr_ROI_bln = strcmp(src_label, ROIs_lbl{ROI_id});
    curr_filters = filters(curr_ROI_bln);
    curr_filters = cat(1,curr_filters{:});
    
    ESI_tmp = time_courses.'*curr_filters.'; % Inverse space transformation
    [ROI_traces(ROI_id,:),S,V] = svds(ESI_tmp,1);
    
end

%% Tell Datalad to allow files to be modified
if exist(derivatives_path, 'dir')
    [status,cmdout] = system('datalad unlock -d '+convertCharsToStrings(BIDS_folder)+' '+ convertCharsToStrings(derivatives_path));
    sprintf(cmdout)
end

%% recreate a Fieldtrip data structure and store it
ROI_data.trial = ROI_traces;
ROI_data.trial = mat2cell(ROI_data.trial, n_ROIs, n_samples*ones(1,n_trials));
ROI_data.elec = ROIs_elec;
ROI_data.trialinfo = data.trialinfo;
ROI_data.label = ROIs_lbl;
ROI_data.fsample = data.fsample;
ROI_data.time = data.time;
ROI_data.sampleinfo = data.sampleinfo;

save(fullfile(derivatives_path,'ROI_data.mat'), 'ROI_data')


%% verify that everything went well
% 
% cfg = [];
% cfg.trials = find(ROI_data.trialinfo == 1);
% avg_FAC = ft_timelockanalysis(cfg, ROI_data);
% cfg = [];
% cfg.trials = find(ROI_data.trialinfo == 0);
% avg_SRC = ft_timelockanalysis(cfg, ROI_data);
% 
% figure;
% cfg = [];
% cfg.layout = ft_prepare_layout(struct('elec', ROI_data.elec));
% ft_multiplotER(cfg,avg_FAC,avg_SRC);
% ft_movieplotER(cfg, avg_FAC)
% 
% cfg = [];
% cfg.channel = 'all';
% cfg.latency = [-0.1 1];
% cfg.method = 'analytic';
% cfg.contrast = [-1 1];
% cfg.statistic = 'indepsamplesT';
% % cfg.numrandomization = 500;
% cfg.correctm         = 'no';
% cfg.alpha            = 0.001;
% cfg.tail             = 0;
% 
% cfg.design = ROI_data.trialinfo+1;
% cfg.ivar = 1;
% 
% stats = ft_timelockstatistics(cfg,ROI_data);
% 

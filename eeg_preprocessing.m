%This script is used to preprocess the EEG data:
%   - import the data 
%   - rereference to A1 to remove the contribute of CMS
%   - downsample to 250 Hz
%   - filter the data (band pass and notch)
%   - segment into trials according to task
%   - (define neihbouring electrodes to be used for interpolation)
%   - identify and remove noisy channels
%   - perform ICA on data
%   - select ICs that contain artefacts 
%   - backproject the data using only the good ICs
%   - substitute the noisy channels by interpolation of the neighbours
%   - average rereference

clc
close all
clear

% make sure you are in the correct folder
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))

addpath(genpath('.')) % add all files and folder to Matlab path

try
    ft_defaults % check if fieldtrip is in Matlab path
catch
    addpath(fullfile('/','home','sinergiasummerschool','Softwares','fieldtrip'))
    ft_defaults
end

%% initialisation
%subj init
BIDS_folder=fullfile('..','..');

task='task-faces';
id_sub=1;
sub_id = sprintf('sub-%02d',id_sub);
filename = fullfile(BIDS_folder, sub_id,'eeg', [sub_id,'_',task,'_eeg.bdf']); %eeg file
event_file=fullfile(BIDS_folder, sub_id,'eeg', [sub_id,'_',task,'_events.tsv']); %event file

%preproc init 
cfg_preproc.resamplefs = 250; % resampling freq
cfg_preproc.hpfreq = 1; % high pass freq
cfg_preproc.bsfreq = [48,52]; %band stop freq / notch freq
cfg_preproc.prestim = 1.5;
cfg_preproc.poststim = 1;
cfg_preproc.rerefmethod = 'avg';



derivatives_path = fullfile(BIDS_folder, 'derivatives', 'eeg_preprocessing', sub_id, 'eeg');
if ~exist(derivatives_path, 'dir')
   mkdir(derivatives_path) 
end


%% preprocess EEG data: reref, downsmpl, filtering
% import the data
cfg.dataset= filename;
cfg.detrend = 'yes';
data = ft_preprocessing(cfg);

%remove STATUS ch and reref to central ch to remove CMS contribute
cfg=struct('channel',{data.label(1:end-1)},'reref','yes','refchannel','A1');
data = ft_preprocessing(cfg,data);

%downsample 
cfg = struct('resamplefs',cfg_preproc.resamplefs);
data_resamp = ft_resampledata(cfg,data);

% filter the data
cfg = struct('hpfilter', 'yes','hpfreq', cfg_preproc.hpfreq, 'bsfilter','yes','bsfreq',cfg_preproc.bsfreq);
data_filt = ft_preprocessing(cfg,data_resamp);

%% read the events.tsc file and downsample the onset sample and duration which were calculated with the original sampling frequency
[event] = ft_read_event(filename);

event(isnan([event(:).sample])) = []; 
for evt = 1:length(event)
    event(evt).sample = round(event(evt).sample/data.fsample*cfg_preproc.resamplefs);
    event(evt).duration = round(event(evt).duration/data.fsample*cfg_preproc.resamplefs);
end


%% segment the data
%identify trials corresponding to faces (1) and scramble (0) 
%-bad epochs are automatically excluded

cfg.hdr.Fs = cfg_preproc.resamplefs;
cfg.hdr.nSamples = data_filt.sampleinfo(2);
cfg.hdr.nTrials = 1;
cfg.event = event;
cfg.tracktimeinfo = 'no';
cfg.trialdef.eventtype='FACES';
cfg.trialdef.eventvalue =1;
cfg.trialdef.prestim = cfg_preproc.prestim; % in seconds
cfg.trialdef.poststim = cfg_preproc.poststim; % in seconds

[cfg_FAC.trl, cfg_FAC.event] = ft_trialfun_general(cfg);

cfg.trialdef.eventtype='SCRAMBLED';
cfg.trialdef.eventvalue =0;
[cfg_SCR.trl, cfg_SCR.event] = ft_trialfun_general(cfg);

%extract trials and exclude trl corresp to 'outliers'
data_seg=segment_trials(cfg_SCR,cfg_FAC,data_filt,event_file,'outliers',1);

%% preprocess EEG data: noisy ch, ICA, avg reref
%prepare neighbours
cfg = struct('channel','all','method','triangulation');
neighbours = ft_prepare_neighbours(cfg, data_seg);

%plot neighbours
cfg=struct('neighbours',neighbours,'elec',data_seg.elec);
ft_neighbourplot(cfg)

%plot layout
layout='biosemi128.lay';
figure;cfg = [];cfg.layout = layout;ft_layoutplot(cfg)

%% read bad channels from *_channels.tsv
chn_filename = strrep(filename,'eeg.bdf','channels.tsv');
chn_tbl = readtable(chn_filename, 'FileType', 'text');
bad_channels = chn_tbl.name(strcmp(chn_tbl.status,'bad'));

% remove bad channels
good_ch=setdiff(data_seg.label,bad_channels);
cfg=[];
cfg.channel=good_ch;
data_no_bad_ch=ft_selectdata(cfg,data_seg);

%store bad channels for json
cfg_preproc.bad_channels=bad_channels;

%% ICA artefact removal

if exist(fullfile(derivatives_path,'ica_component.mat'),'file')
    load(fullfile(derivatives_path,'ica_component.mat'), 'ICs', 'ICs2remove')
else
    % perform ica
    cfg        = [];
    cfg.method = 'runica';
    cfg.channel = {'all' '-A1'}; % compute ICA using all channels except the A1 which is the new reference
    cfg.updatesens = 'no';
    ICs = ft_componentanalysis(cfg, data_no_bad_ch);
    
    %review ICs and select those to remove
    cfg = [];
    cfg.layout = layout; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    cfg.continuous='no';
    cfg.component = 1:20;   % specify the component(s) to plot
    cfg.allowoverlap='yes';
    ft_databrowser(cfg, ICs);
    
    ICs2remove = [1];
    save(fullfile(derivatives_path,'ica_component.mat'),'ICs', 'ICs2remove')
end


%store IC2remove for json
cfg_preproc.ICs2remove=ICs2remove;
cfg=struct('component',ICs2remove);   
cfg.updatesens = 'no';
data_ICA = ft_rejectcomponent(cfg, ICs, data_no_bad_ch);

%% interpolate bad channels and 
cfg = [];
cfg.badchannel=bad_channels;
cfg.method='spline';
cfg.neighbours=neighbours;
data_interp = ft_channelrepair(cfg,data_ICA);

% reset the original order of the channels if changed
orig_labels=data_seg.label;
if ~all(strcmp(data_interp.label, orig_labels))
    [~, idx] = ismember(orig_labels,data_interp.label);
    data_interp.label = data_interp.label(idx);
    data_interp.trial = cellfun(@(trl) trl(idx,:), data_interp.trial, 'uni', 0);
end

%average reference
cfg = struct('reref','yes','refchannel','all','method',cfg_preproc.rerefmethod);
data_reref = ft_preprocessing(cfg, data_interp);

%visualise the final result
cfg=[];
cfg.viewmode='vertical';
cfg.allowoverlap='yes';
ft_databrowser(cfg,data_reref);

%% save to BIDS
clear data
data=data_reref;

%save cleaned data
save(fullfile(derivatives_path,'data.mat'),'data', 'cfg_preproc')

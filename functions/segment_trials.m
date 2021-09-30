function data_seg=segment_trials(cfg_1,cfg_2,data,event_file,event2consider,value2consider)

%This function segments the continuous data according to the values in cfg1
%and cfg2, which excludes the trials with value 'NaN'. Then, the trials
%that have the field 'event2consider' of value 'value2consider' are removed
%from the list as well. Then all the good trials are appended in a single structure

%use as'data_seg=segment_trials(cfg_1,cfg_2,data,event_file,'outliers',1)' if you
%want to remove the trials that are labelled 'outliers=1';

%use as'data_seg=segment_trials(cfg_1,cfg_2,data,event_file,'response_accuracy',0)' if you
%want to remove the trials that are labelled 'response_accuracy=0';


%segment data (NaN trls/bad_epochs excluded)
data_seg_1 = ft_redefinetrial(cfg_1, data);
data_seg_2 = ft_redefinetrial(cfg_2, data);

%read event file
evt_table = readtable(event_file, 'FileType', 'text');

%remove nan index from event_struct
nan_bln = strcmp(evt_table.onset, 'n/a');
evt_table(nan_bln,:) = [];
%find index of trig 2 consider
idx_outliers=find(evt_table.(event2consider) == value2consider); %bad epochs

idx_1=find(evt_table.value==0);
idx_2=find(evt_table.value==1);

%
[bln_1,idx2remove_1]=ismember(idx_outliers,idx_1);
if any(bln_1)
    trl2keep_1=setdiff(1:size(data_seg_1.trialinfo,1),idx2remove_1);
    cfg=struct('trials',trl2keep_1);
    data_seg_redef_1=ft_redefinetrial(cfg,data_seg_1);
else
    data_seg_redef_1=data_seg_1;
end
[bln_2,idx2remove_2]=ismember(idx_outliers,idx_2);
if any(bln_2)
    trl2keep_2=setdiff(1:length(data_seg_2),idx2remove_2);
    cfg=struct('trials',trl2keep_2);
    data_seg_redef_2=ft_redefinetrial(cfg,data_seg_2);
else
    data_seg_redef_2=data_seg_2;
end

%append all the good trials in one single structure
cfg=[];
cfg.keepsampleinfo='yes';
data_seg=ft_appenddata(cfg,data_seg_redef_1,data_seg_redef_2);
function [ICs2remove] = select_ICs2Bremoved(ICs, elec, layout)
% select_ICs2Bremoved open the Fieldtrip mat file containing the
% independent components and asked trhough a GUI which Ics must be removed
%
% Syntax:  [ICs2remove] = select_ICs2Bremoved(filename, elec)
%
% Inputs:
%   filename    - string or cell of strings, fullfile name(s) of the EEG data to
%               preprocess
%   elec        - structure, fieldtrip electrode structure
%   prev_ICs2remove - cell array of scalar, list of ICs (INDICES) that the 
%               user previously chose to remove;
%
% Outputs:
%   ICs2remove  - vector of scalar, with the indices of the
%               components to be removed
%
% Exaxmple:  
%   
%   OR
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Toolbox required: Fieldtrip, A-2_EEG_preprocessing
%
% See also:
% Author: 
%               Nicolas Roehri, Ph.D., neuroscience
%               Isotta Rigoni, MSc
%
% email address:    nicolas.roehri@unige.ch, roehri.nicolas@gmail.com, 
%                   isotta.rigoni@unige.ch.
% Apr 2021; Last revision: 07-Apr-2021

ICs2remove = [];

cfg = [];
cfg.layout = layout; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
cfg.continuous='no';
cfg.ylim        = 1.25 * [-1 1];
cfg.component = 1:20;   % specify the component(s) to plot
cfg.allowoverlap='yes';
ft_databrowser(cfg, ICs);

ft_fig = gcf;
add_rejectICs_button(ft_fig)

set(ft_fig, 'CloseRequestFcn', @cleanup_cb);
  
while ishandle(ft_fig)
    uiwait(ft_fig);
    opt = getappdata(ft_fig, 'opt');
    if opt.cleanup
        cfg = getappdata(ft_fig, 'cfg');
        ICs2remove = cfg.ICs2remove;
        delete(ft_fig);
    end
end

ICs2remove = strrep(ICs2remove, 'runica', '');
ICs2remove = str2double(ICs2remove);

end

function add_rejectICs_button(ft_fig) 
% add Reject components button in ft_databrowser interface
reCompUI = uicontrol('tag', 'labels',  'parent', ft_fig, 'units', 'normalized',...
    'style', 'pushbutton', 'string', 'Reject components', 'Callback',@selection);
reCompUI.Position = [ 0.9100    0.2   0.0800    0.100];
reCompUI.BackgroundColor = [0.5, 1, 0];
end

function selection(h, eventdata)
opt = getappdata(h.Parent, 'opt');
cfg = getappdata(h.Parent, 'cfg');
if ~isfield(cfg, 'ICs2remove')
    cfg.ICs2remove = '';
end

select = match_str(opt.hdr.label, cfg.ICs2remove);
select = select_channel_list(opt.hdr.label, select);
cfg.ICs2remove = opt.hdr.label(select);
setappdata(h.Parent, 'cfg', cfg);
end

function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.cleanup = true;
setappdata(h, 'opt', opt);
uiresume
end % from Fieldtrip
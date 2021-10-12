function marking_win = fepsp_markings(varargin)

% creates gui for user to mark start (green) and end (red) points of a
% response. user can choose to go through all intensities or just one
% intensity per channel (see fast_mark). the number of marker pairs in
% a trace is determined by the number of stimulations derived from the
% stimulus protocol. the output consists of markings for each trace
% separately and markings for the average trace. 
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "GetProtocol.m" for
%                 more info.

% INPUT (optional):
%   intens      - Numeric vector of stimulus intensities used for each
%                 column of traces. Typically provided in units of uA. if
%                 empty will be set as a unit increasing vector.
%   base_path -   String or char. Full path to where the output should be
%                 saved.
%   traces_xlim - Numeric vector of 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty will use protocol default values. See
%                 GetWindowsFromProtocol for more info. 
%                 Default: [].
%   traces_ylim - Numeric vector with 2 elements. voltage limits for
%                 displaying the trace [mV]. these y-limits remain constant
%                 throughout all intensities. If empty will be set
%                 according to the range of traces in highest intensity (excluding stimulus
%                 artifact). 
%                 Default: [].
%   dt          - Non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See GetWindowsFromProtocol for more info.
%                 Default: 2.
%   max_jitter  - Non-negative scalar. Time tolerance for selecting
%                 the response boundaries [ms]. I.e., the precise start
%                 / end points of a response will be the max or min value
%                 within +/- max_jitter of the points selected by the user.
%                 This jitter compensates for small inaccuracies when
%                 marking the peak / trough.
%                 Default: 3.
%   fast_mark   - Logical flag. If true, user will be asked to
%                 mark response boundaries for only the largest intensity.
%                 Those boundaries will be applied to all intensities
%                 without manual inspection. Note that for each intensity,
%                 the precise points of end / start will still be within
%                 +/- max_jitter of the manual selection. 
%                 Default: false.
%
% OUTPUTS:
%   while gui is open the function will return marking_win which is a
%   figure handle to the gui.
%   
%   after the gui window is closed and/or the user clicks "export", the
%   output will consist of the struct markings with the following fields:
% 
%   mark_peaks  - standard cell. 2d cell array of channel (row) x intensity
%                 (column), with each cell containing a 2d numeric mat of
%                 stimulus number (row) x repetition (column). mat values
%                 are the sample of the response peak. note: Any traces
%                 removed by the user will receive nan as their markigns.
%   mark_starts - standard cell. same as mark_peaks but with mat values
%                 representing the sample of the response start.
%   mark_peaks_avg
%               - 2d cell array of channel (row) x intensity (column), with
%                 each cell containing a numeric column vector. vec values
%                 are samples of the response peak, one for each stimulus.
%                 This is similar to "standard cell" (see README) but
%                 inside each cell there is only 1 column for the average
%                 trace instead of a column for every repetition (trace).
%   mark_starts_avg
%               - same as mark_peaks_avg but with vec values representing
%                 the sample of the response start. 
%
% CALL:
%   marking_win = fepsp_markings("traces", <cell array>, "fs", <numeric
%   scalar>, "protocol_id", <string scalar>, "base_path", <folder path>,
%   "intens", <numeric vector>, "traces_xlim", <numeric 2 elements or
%   empty>, "traces_ylim", <numeric 2 elements, or empty>, "dt", <numeric
%   non-negative scalar>, "max_jitter", <numeric non-negative
%   scalar>, "fast_mark", <logical flag>);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Give Option to enter inputs via struct, instead of name-value
p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('traces',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',[],@(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',[],@(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('intens',[],@(x) validateattributes(x,{'numeric'},{'vector'}))
p.addParameter('base_path',pwd,@isfolder)
p.addParameter('traces_xlim',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('traces_ylim',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('dt',2,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))
p.addParameter('max_jitter',3, @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
p.addParameter('fast_mark',false,@(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))

p.parse(varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
intens          = p.Results.intens;
base_path       = p.Results.base_path;
traces_xlim     = sort(p.Results.traces_xlim);
traces_ylim     = sort(p.Results.traces_ylim);
dt              = p.Results.dt;
max_jitter      = p.Results.max_jitter;
fast_mark       = logical(p.Results.fast_mark);

% validate intens
if isempty(intens)
    intens = -(size(traces, 2) : -1 : 1);
elseif numel(intens) ~= size(traces, 2)
    error('number of intensities does not fit size of traces')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stim protocol info and windows
protocol_info = GetProtocol(protocol_id);
windows_info = GetWindowsFromProtocol(protocol_info, fs, dt, traces_xlim);

% initilialize cell arrays of markings
[mark_peaks, mark_starts] =...
    deal(cellfun(@(x) nan(protocol_info.nstim, size(x, 2)), traces, 'UniformOutput', 0));
[mark_peaks_avg, mark_starts_avg] =...
    deal(repmat({nan(protocol_info.nstim, 1)}, size(traces)));

% Convert max_jitter to max index difference from selected time point -
% half in each direction
max_IND_tol = round(((max_jitter/1000)*fs)./2);

% Create avg trace
for iTraces_group = numel(traces):-1:1
    [Chan,Inten] = ind2sub(size(traces),iTraces_group);
    avg_traces(Chan,Inten,:) = mean(traces{iTraces_group},2,'omitnan'); % electrode-by-intensity-by-sample
end

% Create time stamps
Tstamps = -protocol_info.StimLatency:(1/fs)*1000:(protocol_info.TraceRecordingLength-protocol_info.StimLatency);
Tstamps = Tstamps';

% Find the order of intensities
[intens_sorted,intens_order] = sort(intens,'ascend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marking_win = figure();
marking_win.Name = 'fepsp Marking GUI';

% Plot the maximum intensity on channel 1. intens_order(end) is the column
% of the highest intensity
traces_plot = plot(Tstamps,traces{1,intens_order(end)});

% Organize axis limits, from Windows created earlier.
axis tight
xlim(Tstamps(windows_info.DispWindowBase))
if isempty(traces_ylim)
    % Use windows_info.AnalyseWin in order to find the minmax values of the Y
    % axis, without including stim artifact
    yLimit_base = [min(traces{1, intens_order(end)}(windows_info.AnalyseWin,:), [], 'all'),...
        max(traces{1, end}(windows_info.AnalyseWin,:), [], 'all')];
    ylim(yLimit_base*1.1)
else
    % If user defined ylimits, use them
    ylim(traces_ylim)
end

% Label & title nicely :)
xlabel('Time [ms]')
ylabel('Amp [mV]')
title({sprintf('Channel %d @ %duA',1,intens(end))...
    'Move green / red lines to start / peak of each response accordingly'...
    'Press "Enter / return" to save & move to the next intensity'})


% Create Movable Lines.
% Their starting location is based on windows, in order to place them
% "around" each stimulation.
lines_start_pos = Tstamps(windows_info.AnalyseWinBase)';
colors = {'g','r'}; % Prepare for loop - switch between green and red lines
line_rule = {'S','P'}; % Prepare for loop - switch between start (S) and peak (P) lines
% Make lines 20 times bigger then needed to present response. Robust to
% user zoom changes, while maintain label position
lines_Ylim = ylim()'+[-1;1]*(20*max(abs(ylim())));
for iLine = numel(lines_start_pos):-1:1
    rule_IND = mod(iLine,2); % Switch between colours (and S/P marks, accordingly)
    rule_IND(rule_IND == 0) = 2; % Index 0 is bad. Convert it to 2 in order to complete the switch circle
    mark_lines(iLine) = drawline('Position',[repmat(lines_start_pos(iLine),2,1) lines_Ylim],'InteractionsAllowed','translate','color',colors{rule_IND},...
        'Label',sprintf('%d%s t=%.1f',ceil(iLine/2),line_rule{rule_IND},lines_start_pos(iLine)));
    % 'LabelAlpha',0.3 - Not present in 2019a, however recommended in later
    % versions. Helps preventing label from "obscuring" traces. To add it simply
    % add this name value parameter presented here to the "drawline" call
end
Lis = addlistener(mark_lines,'MovingROI',@line_move_manage); %Listener for managing lines movement & label

% Build GUI Buttons & Interactivity:

% Move channel droplist
Ch_strings = split(num2str(1:size(traces,1)));
Ch_button = uicontrol(marking_win,'Style','popupmenu','String',Ch_strings,'Value',1,'Tooltip',sprintf('Channel Selection)'));
% Move intensity droplist
intens_button = uicontrol(marking_win,'Style','popupmenu','String',split(num2str(intens_sorted)),'Value',length(intens_sorted),'Tooltip','Intens Selection');
% Callbacks to move intensity & trace droplists - they are dependent on each
% other, so they need to be defined after they were created
intens_button.Callback = @(~,~) switch_traces_group(marking_win,Ch_button,intens_button,'Int');
Ch_button.Callback = @(~,~) switch_traces_group(marking_win,Ch_button,intens_button,'Ch');
% Export marking to workspace
save_button = uicontrol(marking_win,'Style','pushbutton','String','Export Marks','Tooltip','Export Markings','Callback',@(~,~) export_marks(marking_win));
% Checkbox for fast analysis option
fast_mark_checkbox = uicontrol(marking_win,'Style','checkbox','String','Fast Mark','Value',fast_mark);
% Fix buttons with too small width
save_button.Position(3) = 100;
fast_mark_checkbox.Position(3) = 100;
% Align all button to sit next to each other. If you create any more
% uicontrols, add them here to avoid manual placements
align([Ch_button, intens_button, save_button, fast_mark_checkbox],'Fixed',5,'bottom')

% Create context menu to remove or invert traces from it
cm = uicontextmenu('Callback',@CM_show_nTrace);
uimenu(cm,'Label','TraceNum -') % Using the function trace_trapper, this will mark what plot user right clicked on
uimenu(cm,'Label','Remove Trace','Callback',@(obj,~) remove_trace(obj,marking_win,Ch_button.Value,intens_button.Value));
uimenu(cm,'Label','Invert Trace','Callback',@(obj,~) invert_trace(obj,marking_win,Ch_button.Value,intens_button.Value));
% Note: using 'Callback' property of uicontextmenu is not recommended since
% 2020a. It is recommended to use 'ContextMenuOpeningFcn', its meaning is
% clearer. In order to be compatible to 2019a, we will use
% 'Callback', while remembering this callback run just before uicontextmenu
% is opened.

% Add the context menu to the plots
for iTrace = 1:length(traces_plot)
    traces_plot(iTrace).UIContextMenu = cm;
    traces_plot(iTrace).Tag = num2str(iTrace); % Tag each trace by its column number in the traces group
    traces_plot(iTrace).ButtonDownFcn = {@trace_trapper,cm}; % Improve right click success in finding the right trace
    r = dataTipTextRow('TraceNum',ones(size(traces_plot(iTrace).XData))*iTrace); % Write trace num in datatips
    traces_plot(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
end
% Note: using 'UIContextMenu' property of line is not recommended since
% 2020a. It is recommended to use 'ContextMenu'. 
% In order to be compatible to 2019a, we will use 'UIContextMenu'.

% Add average trace
ax = findobj(marking_win,'type','axes','-depth',1);
hold(ax,'on')
traces_plot(end+1) = plot(Tstamps,squeeze(avg_traces(1,end,:)),'--b','LineWidth',2);
traces_plot(end).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('AvgTrace','');

% Make figure callbacks
marking_win.CloseRequestFcn = @(~,~) close_GUI(marking_win,Lis); % Export marks on window close
marking_win.KeyReleaseFcn =  @(~,evt) next_on_enter(evt,marking_win,Ch_button,intens_button); % Create move trace by pressing "enter"/"return"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save between callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save all data related to the fepsp under 1 struct, to keep it organized
fepsp_data.traces               = traces;
fepsp_data.avg_traces           = avg_traces;
fepsp_data.intens               = intens;
fepsp_data.fs                   = fs;
fepsp_data.protocol_info        = protocol_info;
fepsp_data.Tstamps              = Tstamps;
fepsp_data.windows_info         = windows_info;
fepsp_data.max_IND_tol          = max_IND_tol;
fepsp_data.mark_starts        = mark_starts;
fepsp_data.mark_peaks         = mark_peaks;
fepsp_data.mark_starts_avg    = mark_starts_avg;
fepsp_data.mark_peaks_avg     = mark_peaks_avg;

% Save everything we want to transfer between callbacked under UserData
marking_win.UserData.fepsp_data                     = fepsp_data;
marking_win.UserData.traces_handles.traces_plot     = traces_plot;
marking_win.UserData.traces_handles.cm              = cm;
marking_win.UserData.mark_lines_handles.Lis         = Lis;
marking_win.UserData.mark_lines_handles.mark_lines  = mark_lines;
marking_win.UserData.last_Ch                        = Ch_button.Value;
marking_win.UserData.last_intens                    = intens_button.Value;
marking_win.UserData.fast_mark_checkbox             = fast_mark_checkbox;
marking_win.UserData.intens_order                   = intens_order;
marking_win.UserData.traces_ylim                  = traces_ylim;
marking_win.UserData.base_path                    = base_path;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_markings(marking_win,Ch_button_val,intens_button_val)
% Saving the current markings according to channel & intensity values.
% Before saving, finds the best match within tolerance for each trace.
%
% Inputs:
%   marking_win - Required, figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val -
%                 Required, numeric scalar. Index of the channel to save
%                 to.
%   intens_button_val -
%                 Required, numeric scalar. Index of the intensity to save
%                 to.

% Extract all needed variables form UserData.fepsp_data
fepsp_data       = marking_win.UserData.fepsp_data;
Tstamps          = fepsp_data.Tstamps;
protocol_info    = fepsp_data.protocol_info;
max_IND_tol      = fepsp_data.max_IND_tol;
traces           = fepsp_data.traces;
avg_traces       = fepsp_data.avg_traces;
baseline_window  = fepsp_data.windows_info.BaselineWindow;
intens_order     = marking_win.UserData.intens_order;
analyse_win_base = fepsp_data.windows_info.AnalyseWinBase;

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
traces_col = intens_order(intens_button_val);

% Find the X value of each marking line
lines_pos = [marking_win.UserData.mark_lines_handles.mark_lines.Position];
line_Xpos = lines_pos(1,1:2:end);

% Calculate the closest Point that exist in data
Tstamp_mat = repmat(Tstamps, 1, protocol_info.nstim);
[~,start_IND_base] = min(abs(Tstamp_mat - line_Xpos(1:2:end)), [], 1);
[~,peak_IND_base] = min(abs(Tstamp_mat - line_Xpos(2:2:end)), [], 1);

% Get mean baseline values for each trace & average trace
baseline_traces = mean(traces{Ch_button_val,traces_col}(baseline_window,:),1,'omitnan');
baseline_AVG_trace = mean(avg_traces(Ch_button_val,traces_col,baseline_window),'omitnan');

% Deal with position tolerance (find point in area by max_IND_tol with max
% distance from mean value of baseline)
for iStim = protocol_info.nstim:-1:1
    % Crate tolerance windows for start of the response and for peak:
    tol_window_Start = max(start_IND_base(iStim) - max_IND_tol, analyse_win_base(iStim,1)); % The max makes sure Start tolerance doesn't take stimulus artifact
    tol_window_End = min(start_IND_base(iStim) + max_IND_tol, peak_IND_base(iStim) - max_IND_tol - 1); % The min makes sure Start tolerance doesn't take after peak_IND
    start_tol_window = tol_window_Start:tol_window_End;
    
    tol_window_Start = max(peak_IND_base(iStim) - max_IND_tol, start_IND_base(iStim) + max_IND_tol + 1); % The max makes sure End tolerance doesn't take before start_IND
    tol_window_End = min(peak_IND_base(iStim) + max_IND_tol, analyse_win_base(iStim,2)); % The max makes sure End tolerance doesn't take next stimulus artifact
    peak_tol_window = tol_window_Start:tol_window_End;
    
    % Find the maximal distance from mean baseline in all traces:
    % baseline_traces is a row vector, and inside each traces
    % cell there is a mat in which each column is a different trace
    % (repetition). Therefore, the "minus" operation between the mat & the
    % row vector removes the vector from each row of the matrix, resulting
    % in a distance mat
    distance_mat = abs(traces{Ch_button_val,traces_col}(start_tol_window,:) - baseline_traces);
    [start_values,start_IND(iStim,:)] = max(distance_mat,[],1); 
    distance_mat = abs(traces{Ch_button_val,traces_col}(peak_tol_window,:) - baseline_traces);
    [peak_values,peak_IND(iStim,:)] = max(distance_mat,[],1); 
    % Right now each index refer to the tolerance window.
    % Convert it to full trace
    start_IND(iStim,:) = start_tol_window(1) + start_IND(iStim,:);
    peak_IND(iStim,:) = peak_tol_window(1) + peak_IND(iStim,:);
    % If the minimal value is nan, trace was deleted.
    % However, each IND will mistakenly return 1 + Start of tolerance area.
    % Fix it to NaN.
    start_IND(iStim,isnan(start_values)) = nan;
    peak_IND(iStim,isnan(peak_values)) = nan;
    
    % Find the maximal distance from mean baseline in the average traces:
    distance_vec = abs(avg_traces(Ch_button_val,traces_col,start_tol_window) - baseline_AVG_trace);
    [~,start_IND_AVG(iStim)] = max(distance_vec); 
    distance_vec = abs(avg_traces(Ch_button_val,traces_col,peak_tol_window) - baseline_AVG_trace);
    [~,peak_IND_AVG(iStim)] = max(distance_vec); 
    % Right now each index refer to the tolerance window.
    % Convert it to full trace
    start_IND_AVG(iStim) = start_tol_window(1) + start_IND_AVG(iStim);
    peak_IND_AVG(iStim) = peak_tol_window(1) + peak_IND_AVG(iStim);
   
end

% Save the markings to the right location in each matching
% UserData.fepsp_data variable
marking_win.UserData.fepsp_data.mark_starts{Ch_button_val,traces_col}     = start_IND;
marking_win.UserData.fepsp_data.mark_peaks{Ch_button_val,traces_col}      = peak_IND;
marking_win.UserData.fepsp_data.mark_peaks_avg{Ch_button_val,traces_col}  = peak_IND_AVG';
marking_win.UserData.fepsp_data.mark_starts_avg{Ch_button_val,traces_col} = start_IND_AVG';
end
function switch_traces_group(marking_win,Ch_button,intens_button,Parm)
% Bureaucratic function - takes the window, current channel & current intensity info,
% and what type of transfer user want to make, and switch (replot) to the
% matching traces group. Save markings before moving.
%
% Inputs:
%   marking_win - Required, figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button   - Required, uicontrol handle scalar. The uicontrol in
%                 charge of changing between channels.
%   intens_button - 
%                 Required, uicontrol handle scalar. The uicontrol in
%                 charge of changing between intensities.
%   Parm        - Required, char vector, either 'Next', 'Ch' or 'Int'.
%                 What type of transfer to perform:
%                   Next - move to next traces group (same channel with
%                          lower intensity if available, else next channel
%                          in a rising circular manner - 1 to 2 to 3 back
%                          to 1 - with the highest intensity).
%                   Ch   - Move to user specific channel (according to Ch_button
%                          droplist value), do not change the intensity.
%                   Int  - Move to user specific intensity (according to
%                          intens_button droplist value), do not change the channel.

switch Parm
    case 'Next'
        % Quick move to the next traces group
        
        % Fast analysis case
        if marking_win.UserData.fast_mark_checkbox.Value
            % Save current marks for all intensities
            for iIntens = 1:numel(marking_win.UserData.fepsp_data.intens)
                save_markings(marking_win,Ch_button.Value,iIntens)
            end
            
            % Change intensity to smallest intensity to cause transfer to
            % next channel
            intens_button.Value = 1;
        else
            
            % Save the existing markings before moving traces group
            save_markings(marking_win,Ch_button.Value,intens_button.Value)
        end
        
        % Move to next group
        if intens_button.Value == 1
            % We are on the smallest intensity. Next group requires a channel switch.
            
            % Check if all marks have been already filled - if so, prompt
            % export & close GUI. All mean traces marks are not NaN mean
            % that all marks were filled 
            mark_peaks_avg = [marking_win.UserData.fepsp_data.mark_peaks_avg{:}];            
            if all(~isnan(mark_peaks_avg),'all')
                Status = questdlg('Finished All Channels & Intens. Exit?','Finished Main','Yes','No','Yes');
                switch Status
                    case  "Yes"
                        close_GUI(marking_win,marking_win.UserData.mark_lines_handles.Lis)
                        return
                    case []
                        % Unclear answer - cancel move
                        msgbox('User Cancel, Aborting move')
                        
                        % Flash droplists colours green, to confirm to user
                        % change was cancelled successfully
                        Ch_button.BackgroundColor = 'r';
                        intens_button.BackgroundColor = 'r';
                        pause(0.2)
                        Ch_button.BackgroundColor = 'w';
                        intens_button.BackgroundColor = 'w';
                        return
                end
            end
            
            % Correct droplists
            if Ch_button.Value == size(marking_win.UserData.fepsp_data.traces,1)
                % If this is the last channel, return to the first 
                Ch_button.Value = 1;
            else
                % Move to next channel
                Ch_button.Value = Ch_button.Value+1;
            end
            % Move to the highest intensity
            intens_button.Value = length(intens_button.String);
            
            % Plot the new traces group - done by getting the new channel &
            % intensity from the droplists
            plot_traces_group(marking_win,Ch_button.Value,intens_button.Value)
            
            % Flash droplists colours green, to confirm to user change was
            % done successfully
            Ch_button.BackgroundColor     = 'g';
            intens_button.BackgroundColor = 'g';
            pause(0.2)
            Ch_button.BackgroundColor     = 'w';
            intens_button.BackgroundColor = 'w';
        else
            % We are not on the smallest intensity, move to a smaller one
            intens_button.Value = intens_button.Value-1;
            
            % Plot the new traces group - done by getting the new channel &
            % intensity from the droplists
            plot_traces_group(marking_win,Ch_button.Value,intens_button.Value)
            
            % Flash droplist colour green, to confirm to user change was
            % done successfully
            intens_button.BackgroundColor = 'g';
            pause(0.2)
            intens_button.BackgroundColor = 'w';
        end
    case 'Ch'
        % Change only the channel
        
        % Offer the user to save marking before moving to other trace group
        Status = questdlg('Save current markings?','Save current?','Yes','No','Yes');
        switch Status
            % If 'Yes', save. If 'No', continue without saving. If user
            % cancelled (clicked X), it is unclear answer - cancel move.
            case 'Yes'
                save_markings(marking_win,marking_win.UserData.last_Ch,intens_button.Value)
            case ''
                % Notify user
                msgbox('User Cancel, Aborting move')
                % Restore droplist value to what it was before droplist use
                Ch_button.Value = marking_win.UserData.last_Ch;
                
                % Flash droplist colour red, to confirm to user change was
                % cancelled successfully
                Ch_button.BackgroundColor = 'r';
                pause(0.2)
                Ch_button.BackgroundColor = 'w';
                return
        end
        
        % Plot the new traces group - done by getting the new channel &
        % intensity from the droplists
        plot_traces_group(marking_win,Ch_button.Value,intens_button.Value)
        
        % Sometime this fix "focus stays on droplist & cause enter key to not work" issue
        figure(marking_win);
        Ch_button.Enable = 'off';
        drawnow
        Ch_button.Enable = 'on';
        
    case 'Int'
        % Change only the intensity
        
        % Offer the user to save marking before moving to other trace group
        Status = questdlg('Save current Trace?','Save current?','Yes','No','Yes');
        switch Status
            % If 'Yes', save. If 'No', continue without saving. If user
            % cancelled (clicked X), it is unclear answer - cancel move.
            case 'Yes'
                save_markings(marking_win,Ch_button.Value,marking_win.UserData.last_intens)
            case ''
                % Notify user
                msgbox('User Cancel, Aborting move')
                
                % Restore droplist value to what it was before droplist use
                intens_button.Value = marking_win.UserData.last_intens;
                
                % Flash droplist colour red, to confirm to user change was
                % cancelled successfully
                intens_button.BackgroundColor = 'r';
                pause(0.2)
                intens_button.BackgroundColor = 'w';
                return
        end
        
        % Plot the new traces group - done by getting the new channel &
        % intensity from the droplists
        plot_traces_group(marking_win,Ch_button.Value,intens_button.Value)
        
        % Sometime this fix "focus stays on droplist & cause enter key to not work" issue
        figure(marking_win);
        intens_button.Enable = 'off';
        drawnow
        intens_button.Enable = 'on';
end

% Save current channel & intensity numbers - enable saving marking
% & cancel when changing trace group using droplists, as their
% value is changed before callback
marking_win.UserData.last_Ch     = Ch_button.Value;
marking_win.UserData.last_intens = intens_button.Value;
end
function plot_traces_group(marking_win,wanted_Ch,wanted_inten)
% Change the traces group plotted in marking_win according to given channel
% & intensity indices
%
% Inputs:
%   marking_win - Required, figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   wanted_Ch   - Required, numeric scalar. Index of the channel to plot from.
%   wanted_inten -
%                 Required, numeric scalar. Index of the intensity to plot from.

% Extract all needed variables form UserData.fepsp_data
fepsp_data      = marking_win.UserData.fepsp_data;
Tstamps         = fepsp_data.Tstamps;
traces          = fepsp_data.traces;
windows_info    = fepsp_data.windows_info;
avg_traces      = fepsp_data.avg_traces;
intens          = fepsp_data.intens;
intens_order    = marking_win.UserData.intens_order;
traces_ylim   = marking_win.UserData.traces_ylim;

% wanted_inten is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
traces_col = intens_order(wanted_inten);

% Get the axes
ax = findobj(marking_win,'type','axes','-depth',1);

% Return zoom to default
zoom out

% Remove old traces group
delete(marking_win.UserData.traces_handles.traces_plot)

% Plot the new trace group
traces_plot = plot(ax,Tstamps,traces{wanted_Ch,traces_col});

% Add the context menu to the plots
for iTrace = 1:length(traces_plot)
    traces_plot(iTrace).UIContextMenu = marking_win.UserData.traces_handles.cm;
    traces_plot(iTrace).Tag = num2str(iTrace); % Tag each trace by its column number in the trace cell
    traces_plot(iTrace).ButtonDownFcn = {@trace_trapper,marking_win.UserData.traces_handles.cm}; % Improve right click success in finding the right trace
    r = dataTipTextRow('TraceNum',ones(size(traces_plot(iTrace).XData))*iTrace); % Write trace num in datatips
    traces_plot(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
end

% Add average trace
traces_plot(end+1) = plot(Tstamps,squeeze(avg_traces(wanted_Ch,traces_col,:)),'--b','LineWidth',2);
traces_plot(end).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('AvgTrace','');

% Change title to match the new traces group
title(ax,{sprintf('Channel %d @ %duA',wanted_Ch,intens(traces_col))...
    'Move green / red lines to start / peak of each response accordingly'...
    'Press "Enter / return" to save & move to the next intensity'})

% Fix axis limits if necessary:
xlim(Tstamps(windows_info.DispWindowBase))
% Get the needed ylimits to present the current traces group (while ignoring the stim artifact)
traces_group_Ylims = [min(traces{wanted_Ch, traces_col}(windows_info.AnalyseWin,:), [], 'all'),...
    max(traces{wanted_Ch, traces_col}(windows_info.AnalyseWin,:), [], 'all')]*1.1;

if ~isempty(traces_ylim)
    % If user defined ylimits, keep them
    ylim(traces_ylim)
else
    if marking_win.UserData.last_Ch~= wanted_Ch
        % Channel was switch - inherited ylimits are irrelevant. Compare
        % ylimits of biggest intensity to current traces group- take the
        % biggest of the two
        
        % Get the ylimits matching the biggest intensity (while ignoring the
        % stim artifact). intens_order(end) is the column of the highest intensity
        high_intens_Ylims = [min(traces{wanted_Ch, intens_order(end)}(windows_info.AnalyseWin,:), [], 'all'),...
            max(traces{wanted_Ch, intens_order(end)}(windows_info.AnalyseWin,:), [], 'all')]*1.1;
        % Compare and take the biggest
        ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
    else
        % Compare ylimit inherited from last traces group to
        % current - take the biggest of the two.
        
        % Get the ylimits inherited
        inherited_Ylims = ylim();
        % Compare and take the biggest
        ylim([min(traces_group_Ylims(1),inherited_Ylims(1)) max(traces_group_Ylims(2),inherited_Ylims(2))])
    end
end

% Make lines 20 times bigger then current ylimits, by adding. Robust to
% user zoom changes, while maintain label position
lines_Ylim = ylim()'+[-1;1]*(20*max(abs(ylim())));
% Change markers lines to match new ylimits
for aa = 1:length(marking_win.UserData.mark_lines_handles.mark_lines)
    marking_win.UserData.mark_lines_handles.mark_lines(aa).Position(:,2) = lines_Ylim;
end

% Save the handle to the new plots
marking_win.UserData.traces_handles.traces_plot = traces_plot;
end
function remove_trace(obj,marking_win,Ch_button_val,intens_button_val)
% Remove selected trace, and all data connected to it. It uses the context
% menu to find which trace, the window for UserData, and current channel &
% intensity indices in order to remove the data from the correct location
% and redraw the traces group (with the bad trace removed)
%
% Inputs:
%   obj         - Required, uimenu scalar. The caller object for
%                 this callback function.
%   marking_win - Required, figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val - 
%                 Required, numeric scalar. Index of the channel to plot from.
%   intens_button_val - 
%                 Required, numeric scalar. Index of the intensity to plot from.

% Get the trace User want to remove, from its plotted line. It is saved
% under the uicontextmenu UserData, the parent of the uimenu that called this
% callback
Trace = obj.Parent.UserData;
if ~isa(Trace,'matlab.graphics.chart.primitive.Line')
    % If no trace was found, give nice error to user, and cancel
    errordlg({'For some reason, no trace have been capture.' 'Please Try again'})
    return
end
% Get the trace number (its column number in the trace cell) from its tag
nTrace = str2double(Trace.Tag);

% Make sure user want to remove trace - it is not a mistake
Status = questdlg({'Are you sure you want to remove trace?' 'This cannot be undone, and will effect the output struct'},'Removal Confirm','Yes','No','No');
if ismember(Status,{'No',''})
    return
end

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
intens_order = marking_win.UserData.intens_order;
traces_col = intens_order(intens_button_val);

% Delete the trace data
marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace)        = nan;
marking_win.UserData.fepsp_data.mark_peaks{Ch_button_val,traces_col}(:,nTrace)  = nan;
marking_win.UserData.fepsp_data.mark_starts{Ch_button_val,traces_col}(:,nTrace) = nan;
% Recalculate mean trace, for changing display later
marking_win.UserData.fepsp_data.avg_traces(Ch_button_val,traces_col,:)            = mean(marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col},2,'omitnan');

% Redraw the traces - the deleted trace won't be present.
plot_traces_group(marking_win,Ch_button_val,intens_button_val)

if isempty(marking_win.UserData.traces_ylim)
    % If user defined ylimits, they are fixed, and will not be changed. Else:
    % Try to shrink the ylimits, however no smaller than this channel max
    % intens ylimits or the needed to view the traces
    traces = marking_win.UserData.fepsp_data.traces;
    responces_win =  marking_win.UserData.fepsp_data.windows_info.AnalyseWin;
    high_intens_Ylims = [min(traces{Ch_button_val, intens_order(end)}(responces_win,:), [], 'all'),...
        max(traces{Ch_button_val, intens_order(end)}(responces_win,:), [], 'all')]*1.1;
    traces_group_Ylims = [min(traces{Ch_button_val,traces_col}(responces_win,:), [], 'all'),...
        max(traces{Ch_button_val,traces_col}(responces_win,:), [], 'all')]*1.1;
    % Compare and take the biggest
    ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
end

end
function invert_trace(obj,marking_win,Ch_button_val,intens_button_val)
% Invert a selected trace. It uses the context menu to find which trace, the
% window for UserData, and current channel & intensity indices in order to
% invert the data from the correct location and redraw the traces group
%
% Inputs:
%   obj         - Required, uimenu scalar. The caller object for
%                 this callback function.
%   marking_win - Required, figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val - 
%                 Required, numeric scalar. Index of the channel to plot from.
%   intens_button_val - 
%                 Required, numeric scalar. Index of the intensity to plot from.

% Get the trace User want to invert, from its plotted line. It is saved
% under the uicontextmenu UserData, the parent of the uimenu that called this
% callback
trace_line_obj = obj.Parent.UserData;
if ~isa(trace_line_obj,'matlab.graphics.chart.primitive.Line')
    errordlg({'For some reason, no trace have been capture.' 'Please Try again'})
    return
end
nTrace = str2double(trace_line_obj.Tag);

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
intens_order = marking_win.UserData.intens_order;
traces_col = intens_order(intens_button_val);

% Invert the trace data
marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace) = -marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace);

% Recalculate average trace
marking_win.UserData.fepsp_data.avg_traces(Ch_button_val,traces_col,:) = mean(marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col},2,'omitnan');

% Redraw trace group
plot_traces_group(marking_win,Ch_button_val,intens_button_val)

if isempty(marking_win.UserData.traces_ylim)
    % If user defined ylimits, they are fixed, and will not be changed. Else:
    % Try to shrink the ylimits, however no smaller than this channel max
    % intens ylimits or the needed to view the traces
    traces = marking_win.UserData.fepsp_data.traces;
    responces_win =  marking_win.UserData.fepsp_data.windows_info.AnalyseWin;
    high_intens_Ylims = [min(traces{Ch_button_val, intens_order(end)}(responces_win,:), [], 'all'),...
        max(traces{Ch_button_val, intens_order(end)}(responces_win,:), [], 'all')]*1.1;
    traces_group_Ylims = [min(traces{Ch_button_val,traces_col}(responces_win,:), [], 'all'),...
        max(traces{Ch_button_val,traces_col}(responces_win,:), [], 'all')]*1.1;
    % Compare and take the biggest
    ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
end

end
function export_marks(marking_win)
% Export Marking from GUI UserData to base workspace
%
% Inputs:
%   marking_win      - Required, figure scalar, window of fepsp_markings
%

% Pack All Outputs in 1 struct, for easy export
fepsp_markings.mark_peaks = marking_win.UserData.fepsp_data.mark_peaks;
fepsp_markings.mark_starts = marking_win.UserData.fepsp_data.mark_starts;
fepsp_markings.mark_peaks_avg = marking_win.UserData.fepsp_data.mark_peaks_avg;
fepsp_markings.mark_starts_avg = marking_win.UserData.fepsp_data.mark_starts_avg;

% Save to file located in base_path:

% Create file name from base_path
base_path = marking_win.UserData.base_path;
[~,base_name] = fileparts(base_path);
marking_file = [base_path filesep base_name '_fepsp_markings.mat'];
if exist(marking_file,'file')
    % Normal file name is taken. Ask user what to do
    answer = questdlg(sprintf('file "%s" already exist. Do you wish to overwrite it, or a save new file with a numeric suffix "(x)"?\n{Cancelling will save a copy}',marking_file),...
        'Overwrite file?','Overwrite','Create Copy','Create Copy');
    if isempty(answer) || strcmp(answer,'Create Copy')
        % User did not want to overwrite. Create a file name with a numeric
        % suffix that won't overwrite anything
        
        % Initialize
        numeric_sufx = 1;
        marking_file = insertBefore(marking_file,'.mat',sprintf('(%d)',numeric_sufx));
        % Replace suffix until we will find a good number
        while exist(marking_file,'file')
            numeric_sufx = numeric_sufx + 1;
            marking_file = replace(marking_file,sprintf('(%d)',numeric_sufx-1),sprintf('(%d)',numeric_sufx));
        end
    end
end

% Save - tell user where file was saved
fprintf('Saving markings as "%s"\n',marking_file)
save(marking_file,'-struct','fepsp_markings')

end
function close_GUI(marking_win,Lis)
% Simple wrapper callback:
% Save everything before closing GUI. Delete listeners - they should be
% removed when lineROI are deleted, but deleting them explicitly helps in
% some bug cases.
%
% Inputs:
%   marking_win      - Required, figure scalar, window of fepsp_markings
%   Lis         - Required, listeners vector for lineROI event "MovingROI"

export_marks(marking_win)
delete(Lis)
delete(marking_win)

end
function line_move_manage(obj,evt)
% Function for managing lines movement.
% Change the time written in lineROI label when moving it.
% Prevent moving the lines along the Y axis.
%
% Input
%   obj         - Required, lineROI object that generated the event
%                 "MovingROI".
%   evt         - Required, "MovingROI" event, generated in response to
%                 user moving the lineROI.

% Get the new position of the lineROI
XPos = evt.CurrentPosition(1);
% Change the label to include it, after the '=' sign
obj.Label = [obj.Label(1:(find(obj.Label == '='))), sprintf('%.1f',XPos)];

% Prevent moving on the Y axis, by restoring previous position
obj.Position(:,2) = evt.PreviousPosition(:,2);
end
function next_on_enter(evt,marking_win,Ch_button,intens_button)
% Cause movement to the next traces group when figure window is in focus,
% and user pressed "enter"/"return" key.
%
% Inputs:
%   evt         - Required, the keypress event generated by the figure.
%   marking_win - Required, figure scalar, window of fepsp_markings.
%   Ch_button   - Required, uicontrol handle scalar. The uicontrol in
%                 charge of changing between channels.
%   intens_button - 
%                 Required, uicontrol handle scalar. The uicontrol in
%                 charge of changing between intensities.

% Make sure the key that was pressed was "enter"/"return" - its number is 13
if evt.Character == 13
    % If it was, just call switch_traces_group with parameter 'Next'
    switch_traces_group(marking_win,Ch_button,intens_button,'Next')
end
end
function trace_trapper(obj,~,cm)
% Will help the uicontext miss less often, in the cost of hover datatip
% due to some unknown reason (datacrusermode will still work as usual however)
%
% Input
%   obj         - line obj scalar, of the trace that was clicked
%   ~           - The click event generated by the clicked trace
%   cm          - uicontextmenu scalar, opened when user right-click the
%                 trace

% Simply place the trace handle under the context menu UserData for easy
% and reliable access
cm.UserData = obj;
end
function CM_show_nTrace(cm,~)
% Present in trace number 1 of the uimenues what trace the user had right
% clicked on, when opening the uicontextmenu.
%
%   cm          - uicontextmenu scalar, opened when user right-click the
%                 trace, object calling this callback
%   ~           - The generated event from opening the uicontextmenu

% Extract the tag of the trace from the uicontext menu UserData (it was
% saved there by trace_trapper) and switch in the label of the uimenu what after
% the '-' character
cm.Children(3).Label = [extractBefore(cm.Children(3).Label,'-') '- ' num2str(cm.UserData.Tag)];
end
function marking_win = fepsp_markings(varargin)
% Create GUI in order to extract marking of the start and end of
% response for each trace. Returns handle to the window of the GUI when
% function finish function run, and upon closing the GUI saves a file
% with all the markings.
%
% INPUTS:
%   traces      - 2d cell array of channel (row) x intensity (column) with
%                 each cell containing a 2d numeric mat of sample (row) x
%                 repetition (trace, column). This is similar to "Package
%                 Standard Cell Organization" (see also README), but inside
%                 each cell, rows are samples instead of stimulus number.
%                 This is the organized raw data.
%   fs          - Required. Numeric scalar. Sampling frequency [Hz].
%   protocol_id - Required. String scalar or char vector. ID of stimulation
%                 protocol used, for example "io","stp" or "custom" (for
%                 manually writing a new protocol). 
%                 See "GetProtocol" for more info.
%   intens      - Optional. Numeric vector describing the stimulus
%                 intensity used for each column of traces.
%                 Unless traces_ylimits are entered by user, they are
%                 determined automatically for each channel, based on the
%                 traces group with the highest (NOT absolute) intensity.
%                 Default: An increasing vector, from minus number of
%                 columns in traces up to minus 1. For example if traces
%                 has 3 columns, the default intens will be the vector:
%                 [-3,-2,-1].
%   base_folder - Optional. The folder in which the whole project is
%                 managed. The name of this folder will be present in every
%                 file of this project, and if any required input is
%                 missing, we will try to load it from file named in the
%                 base_folder. The file must be named "<base_folder_name>_fepsp_traces.mat"
%                 (for example if base folder path is "C:\Users\Slutsky", the
%                 file has to be named "Slutsky_fepsp_traces.mat").
%                 The file must contain the required inputs arguments
%                 ("traces","fs" & "protocol_id").
%                 See README for more information.
%                 Default: MATLAB's current folder, at the time the function
%                 was called. I.e., if when the function was called the
%                 MATLAB was on "Folder1", and before closing GUI the folder
%                 was changed to "Folder2" - base_folder will still be
%                 "Folder1" in the default case, and output will save there.
%   traces_Xlimit - 
%                 Optional. Numeric vector with 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty ( [] ) will use protocol default values. 
%                 See GetWindowsFromProtocol for more info.
%                 Default: Empty ( [] ).
%   traces_Ylimit -
%                 Optional. Numeric vector with 2 elements. Constant voltage
%                 limits for displaying the trace, in [mV]. Those ylimits
%                 are constant - they remine the default values when
%                 moving between traces groups.
%                 If empty ( [] ) will choose ylimits automatically: Take the
%                 limits needed to present the traces group in current
%                 channel with the highest intensity (excluding stimulus
%                 artifact). Limit will expend if needed for a specific
%                 traces group.
%                 For example, if intens is the vector [20,-100,90,3],
%                 ylimits will be taken from column 3 for each channel.
%                 Automat ylimits use a "window" to exclude stimulus
%                 artifact. See GetWindowsFromProtocol for more info.
%                 Default: Empty ( [] ).
%   dt          - Optional. Non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See GetWindowsFromProtocol for more info.
%                 Default: 2.
%   max_time_tol -
%                 Optional. Non-negative scalar. Time tolerance for selecting
%                 the response boundaries [ms]. I.e., the precise start
%                 / end points of a response will be the max or min value
%                 within +/- max_time_tol of the points selected by the user.
%                 This â€œjitterâ€? compensates for small inaccuracies when
%                 marking the peak / start, and for between traces
%                 variation in the same traces group.
%                 Default: 3.
%   fast_mark -   Optional. Logical flag. If true, user will be asked to
%                 mark response boundaries for only one (the biggest)
%                 intensity. Those boundaries will be applied to all
%                 intensities without manual inspection. 
%                 Note that for each intensity, the precise points of end /
%                 start will still be within +/- max_time_tol of the manual
%                 selection.
%                 Default: false.
%
% OUTPUTS:
% Upon function finish running:
%   marking_win - figure handle, Window Of the GUI. All function data can
%                 be found under UserData.
% After GUI window was closed, or when user clicked "export":
%   A mat file that named "<base_folder_name>_fepsp_markings", will be
%   saved. For example if base path is "C:\Users\Slutsky", the saved file
%   will be "Slutsky_fepsp_markings". Inside this file, the following
%   variables will be:
%       marked_peaks -
%                 2d cell array of channel (row) x intensity (column), with
%                 each cell containing a 2d numeric mat of stimulus number
%                 (row) x repetition (trace, column). This is the "Package
%                 Standard Cell Organization" (see also README).
%                 Each value is the sample number of the response peak.
%                 Note: Any traces removed during marking, their marked
%                 peaks will be NaN for all stimulus numbers (rows).
%       marked_starts - 
%                 "Package Standard Cell Organization" (see 'marked_peaks' above).
%                 Each value is the sample number of the response start.
%                 Note: Any traces removed during marking, their marked
%                 starts will be NaN for all stimulus numbers (rows).
%       Avg_marked_peaks -
%                 2d cell array of channel (row) x intensity (column), with
%                 each cell containing a numeric column vector, rows are
%                 stimulus number. This is similar to "Package Standard
%                 Cell Organization" (see also README), but inside each
%                 cell, there is only 1 column - for the average trace -
%                 instead of a column for every repetition (trace). Each
%                 value is the sample number of the response peak, for the
%                 average trace of the matching traces group in "traces".
%       Avg_marked_starts -
%                 Organized the same way as "Avg_marked_peaks" (see above).
%                 Each value is the sample number of the response starts,
%                 for the average trace of the matching traces group in
%                 "traces".
%
% EXAMPLES:
%   1) * Full Function Call *
%   marking_win = fepsp_markings("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>,...
%   "base_folder",<folder path>,"intens",<numeric vector>,"traces_Xlimit",<numeric 2 elements, or empty>,"traces_Ylimit",<numeric 2 elements, or empty>,...
%   "dt",<numeric non-negative scalar>,"max_time_tol",<numeric non-negative scalar>,"fast_mark",<logical flag>)
%   % Filled example:
%   % This example shows function call, while explicitly defining any
%   % optional input default value. Notice that required inputs -
%   % "traces","fs" & "protocol_id" - do not have default value, and an input
%   % example is shown here instead. 
%   marking_win = ...
%   fepsp_markings("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom",...
%   "base_folder",pwd,"intens",[-2,-1],"traces_Xlimit",[],"traces_Ylimit",[],...
%   "dt",2,"max_time_tol",3,"fast_mark",false)
%
%   2) * Defining only required inputs, using optional defaults by omitting *
%   % Any optional input can be set to default value by simply omitting it.
%   % You can omit all or only some of the optional inputs.
%   marking_win = fepsp_markings("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>)
%   % Filled example:
%   % The filled function call in example 1 is the same as simply:
%   marking_win = ...
%   fepsp_markings("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom")
%
%   3) * Loading missing required inputs from files in "base_folder" *
%   % If any required inputs ("traces","fs","protocol_id') are missing, and
%   % base_folder contain a file containing them (as specified in
%   % "base_folder" input description), they will be loaded. 
%   % For example the following code will produce the same results as
%   % the filled example 2:
%   traces = {ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)};
%   fs = 100;
%   protocol_id = "custom";
%   save("C:\User\Slutsky\Slutsky_fepsp_traces.mat","traces","fs","protocol_id")
%   marking_win = fepsp_markings("base_folder","C:\User\Slutsky")
%   % #Note1: You can use optional inputs by specifying them in this case, as always.
%   % #Note2: calling fepsp_markings without any inputs (fepsp_markings() )
%   % will use the current MATLAB folder as is the default for "base_folder".
%
%   4) * Fast Mark*
%   % In order to perform a fast-marking process, simply set "fast_mark"
%   % flag to "true":
%   fepsp_markings(___,"fast_mark",true)
%   % You can do that with any other syntax.
%   % Filled example: (simillar to example 2)
%   marking_win = ...
%   fepsp_markings("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom","fast_mark",true)
%   % You can also use the "Fast Mark" checkbox in the GUI in order to
%   % preform / cancel fast mark after function call
%
% TROUBLESHOOTING
%   I'm calling fepsp_markings from function or a script, and unable to get
%   the markings:
%       Notice that when function finish running, and only return the GUI
%       window to workspace. The markings are available only after GUI was
%       closed (or user clicked export), and only as a saved file. In order
%       to get the marking into your workspace while running the script /
%       function, have your code wait for the GUI to be deleted, and then
%       load the saved file. Note that if file already exist in base_path &
%       user choose to save a copy, saved file name may be different from
%       expected.
%       Code example:
%       
%       waitfor(fepsp_markings("base_folder","C:/User/Folder1"))
%       Markings = load(fullfile("C:/User/Folder1",['Folder1' '_fepsp_markings']));
%       
%   See also fepsp_org2traces, fepsp_analyse


% Note1: is checking if channel changed for ylimits in PlotANew really
%        nessasry? For changing channel when on low intens. Sound like
%        still need to compare to current plot, in case of out of Y axis.
%        See Note 2. Needed. !!!!Done!!!!
% Note2: Ylimit deciding method is bad. choose rules for basis, expanding &
%        shrinking.
%        Option 1: We always need to check for expanding, relative to
%                  inherited (existing) ylim. We only need to shrink if we
%                  inverted or remove, and only to minimal max intn of
%                  this channel. Problem case: we needed to expand, then
%                  in lower inten we inverted, as a resualt shrink back to
%                  smaller than previosly expended.
%        Option 2: We always first try to expand, then try to shrink back
%                  to minimal max intn of this channel. similar Problem
%                  case as option 1, but now not only when inverting.
%        Option 3: Fix ylim to widest of all this channel traces. Need to
%                  refit when inverting or removing. Problem: bad trace in
%                  1 intens will cause problem in all of them, but won't be
%                  very visiable.
%        Option 4: Try to expend every time I move traces group, realtive to inherited Ylim.
%                  Save ylim each time I move traces group (after expanding), then
%                  when invert/remove try to shrink to those. What to do
%                  when moving not via next? Same as moving by next,
%                  expanding / max intens in channel. What to do if we are
%                  on the max intens in channel? same as always, try to
%                  shrink to max channel = try to shrink to saved. Solve
%                  problem from option 1, as I will return to inherited &
%                  not to max intens. Problems: moving after inverting
%                  "locks in" the ylimits, pre-inverted trace is
%                  problematic.
%        implamented option 2. Working, not ideal, but will do for now.
%        !!!!Done!!!!
% Note3: Remove trace renumber traces. This make it hard to find them
%        later, if needed to re-remove for example. Solution: instead of
%        totaly removing data, NaN it. Implamented, check it, Mainly effect
%        on saving & Ylim. !!!!Done!!!
% Note4: Output as struct is not very comfterbale. Maybe check for every
%        var, and jump questdlg if overwrite or save as ans? what about
%        multi ans case? What about calling from function and not cmd? Save
%        and user can load. !!!!DONE!!!!
% Note5: Should Ylimit take into accuont baseline as well? problem overflow
%        to lineROI initial pos & GetWindowsFromProtocol func !!!!Nope!!!!
% Note6: Itruduce tolerance to start marks. Do not choose maxima/minima for
%        start / end marking. Options:
%        1) Search for max value in tolerance windows under absolute - does
%           not help, everything may be posative / negative. 
%        2) Calculate prominences for the whole anaylse window, both normal
%           & inverted. Take the peak with maximal prominance (from both
%           cases), and try to get the point closest to it in tolerance
%           window.Problem - noisy data, some peaks my have lower
%           prominance then they should as they are "riding" other extrima.
%        3) For each extrima in 1 response, find absolote distance from
%           mean baseline. Find the point with the maximum distance. Take
%           the point closest to it in tolerance window.
%        4) For each point in tolerance window, find absolote distance from
%           mean baseline. Find the point with the maximum distance.
%           Seriosly, how didn't I think about this ages ago?!
%        Implamented (4). Check goodness of fit. fit a littel worst the
%        max / min, but doable. !!!!Done!!!!
% Note7: What will happen if traces are not the same length between cells?
%        Tstamp won't match them. Validate length of each trace match
%        Tstamps: (TraceRecordingLength./(1000/fs)) + 1
% Note8: I'm not sure "End of the response" & "Start of the response" is a
%        good nomenclature. Changed to "Peak" & "Start". Speak about trough
%        in README.
% Note9: Merge AVG markings into all data, as end+1 trace, at least for
%        output. Pro: reduce number of inputs to each function. Con: Make
%        calculation on anything more confusing to user (always remember to
%        remove last trace). !!!!Nope!!!!
% Note10: Should we prevent points to be outside of AnalyseWin? YES,
%         !!!!Done!!!!
% Note11: Weird results case: Saving traces group markers, then returning
%         to it. Removing or inverting trace does change mean trace, but do
%         not currect marking for it. Add warning to "save markings"
%         prompt? Or is it clear that you need to save when you make a
%         change?
% Note12: Give option for tolarance method - far from baseline, start-min
%         peak-max, start-max peak-min. Do it by seperating Start & Peak
%         refernce baseline - mean baseline for far from baseline, max
%         trace value for min & min trace value for max. Always search for
%         the absolute diffrence between points in tolarance & baseline.
% Note13: Fix start_tol_window is empty due to too big dt. Prevent lines
%         getting outside analysis window?
% Note14: assuming start is on the edge of the responce - close to the
%         baseline. It'll always try to run away?

%% Input Parser - Define & validate All Inputs

%Give Option to enter inputs via struct, instead of name-value
p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

% Data arguments
p.addParameter('traces',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',[],@(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',[],@(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('intens',[],@(x) validateattributes(x,{'numeric'},{'vector'}))

% Function management arguments
p.addParameter('base_folder',pwd,@isfolder)
p.addParameter('traces_Xlimit',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('traces_Ylimit',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('dt',2,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))
p.addParameter('max_time_tol',3, @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
p.addParameter('fast_mark',false,@(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))

p.parse(varargin{:})

% Data inputs
traces        = p.Results.traces;
fs            = p.Results.fs;
protocol_id   = p.Results.protocol_id;
intens        = p.Results.intens;

% Function management arguments
base_folder   = p.Results.base_folder;
traces_Xlimit = sort(p.Results.traces_Xlimit);
traces_Ylimit = sort(p.Results.traces_Ylimit);
dt            = p.Results.dt;
max_time_tol  = p.Results.max_time_tol;
fast_mark     = logical(p.Results.fast_mark);

% Make sure all required inputs were given
req_inputs = {'traces','fs','protocol_id'};
req_logIND = ismember(req_inputs,p.UsingDefaults);
if any(req_logIND)
    % If there are missing inputs, try to load them from '_fepsp_traces'
    [~,base_name] = fileparts(base_folder);
    traces_file = [base_folder filesep base_name '_fepsp_traces.mat'];
    if exist(traces_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs{req_logIND}),traces_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(traces_file,req_inputs{req_logIND})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(traces,{'cell'},{'2d'},'fepsp_markings','traces')
        validateattributes(fs,{'numeric'},{'scalar'},'fepsp_markings','fs')
        validateattributes(protocol_id,{'string','char'},{'scalartext'},'fepsp_markings','protocol_id')
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''traces'',''fs'',''protocol_id''} are required.'],sprintf('%s, ',req_inputs{req_logIND}),traces_file)
    end
end


% Validate intens
if isempty(intens)
    % If intens were not given, create mock intens in the negative domain
    intens = -(size(traces,2):-1:1);
elseif numel(intens) ~= size(traces,2)
    error('Number of inten does not match number of intensities in traces (number of columns in the cell array)')
end
%% Crete helper variables

% Get protocol struct info from protocol_id
protocol_info = GetProtocol(protocol_id);

% Windows mark important area in the trace - baseline (before first
% stimulus), display (xlimits) and response area (expected response area,
% for excluding stimulus artifacts).
% See GetWindowsFromProtocol for more info
windows_info = GetWindowsFromProtocol(protocol_info,fs,dt,traces_Xlimit);

% Preallocate Marked Times outputs from GUI - mainly to define their size
% marked_peaks,marked_starts - 
% cell array channel-by-intensity, inside double stim-by-trace
[marked_peaks,marked_starts] = deal(cellfun(@(x) nan(protocol_info.nstim,size(x,2)),traces,'UniformOutput',0));
% Avg_marked_peaks,Avg_marked_starts - 
% cell array channel-by-intensity, inside double stim (rows) col vec
[Avg_marked_peaks,Avg_marked_starts] = deal(repmat({nan(protocol_info.nstim,1)},size(traces)));

% Convert max_time_tol to max index difference from selected time point -
% half in each direction
max_IND_tol = round(((max_time_tol/1000)*fs)./2);

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
%% Build GUI
marking_win = figure();
marking_win.Name = 'fepsp Marking GUI';

% Plot the maximum intensity on channel 1. intens_order(end) is the column
% of the highest intensity
traces_plot = plot(Tstamps,traces{1,intens_order(end)});

% Organize axis limits, from Windows created earlier.
axis tight
xlim(Tstamps(windows_info.DispWindowBase))
if isempty(traces_Ylimit)
    % Use windows_info.AnalyseWin in order to find the minmax values of the Y
    % axis, without including stim artifact
    yLimit_base = [min(traces{1, intens_order(end)}(windows_info.AnalyseWin,:), [], 'all'),...
        max(traces{1, end}(windows_info.AnalyseWin,:), [], 'all')];
    ylim(yLimit_base*1.1)
else
    % If user defined ylimits, use them
    ylim(traces_Ylimit)
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

%% Save data to transfer between callbacks

% Save all data related to the fepsp under 1 struct, to keep it organized
fepsp_data.traces               = traces;
fepsp_data.avg_traces           = avg_traces;
fepsp_data.intens               = intens;
fepsp_data.fs                   = fs;
fepsp_data.protocol_info        = protocol_info;
fepsp_data.Tstamps              = Tstamps;
fepsp_data.windows_info         = windows_info;
fepsp_data.max_IND_tol          = max_IND_tol;
fepsp_data.marked_starts        = marked_starts;
fepsp_data.marked_peaks         = marked_peaks;
fepsp_data.Avg_marked_starts    = Avg_marked_starts;
fepsp_data.Avg_marked_peaks     = Avg_marked_peaks;

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
marking_win.UserData.traces_Ylimit                  = traces_Ylimit;
marking_win.UserData.base_folder                    = base_folder;
end

%% Helper Functions
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
marking_win.UserData.fepsp_data.marked_starts{Ch_button_val,traces_col}     = start_IND;
marking_win.UserData.fepsp_data.marked_peaks{Ch_button_val,traces_col}      = peak_IND;
marking_win.UserData.fepsp_data.Avg_marked_peaks{Ch_button_val,traces_col}  = peak_IND_AVG';
marking_win.UserData.fepsp_data.Avg_marked_starts{Ch_button_val,traces_col} = start_IND_AVG';
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
            Avg_marked_peaks = [marking_win.UserData.fepsp_data.Avg_marked_peaks{:}];            
            if all(~isnan(Avg_marked_peaks),'all')
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
traces_Ylimit   = marking_win.UserData.traces_Ylimit;

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

if ~isempty(traces_Ylimit)
    % If user defined ylimits, keep them
    ylim(traces_Ylimit)
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
marking_win.UserData.fepsp_data.marked_peaks{Ch_button_val,traces_col}(:,nTrace)  = nan;
marking_win.UserData.fepsp_data.marked_starts{Ch_button_val,traces_col}(:,nTrace) = nan;
% Recalculate mean trace, for changing display later
marking_win.UserData.fepsp_data.avg_traces(Ch_button_val,traces_col,:)            = mean(marking_win.UserData.fepsp_data.traces{Ch_button_val,traces_col},2,'omitnan');

% Redraw the traces - the deleted trace won't be present.
plot_traces_group(marking_win,Ch_button_val,intens_button_val)

if isempty(marking_win.UserData.traces_Ylimit)
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

if isempty(marking_win.UserData.traces_Ylimit)
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
fepsp_markings.marked_peaks = marking_win.UserData.fepsp_data.marked_peaks;
fepsp_markings.marked_starts = marking_win.UserData.fepsp_data.marked_starts;
fepsp_markings.Avg_marked_peaks = marking_win.UserData.fepsp_data.Avg_marked_peaks;
fepsp_markings.Avg_marked_starts = marking_win.UserData.fepsp_data.Avg_marked_starts;

% Save to file located in base_folder:

% Create file name from base_folder
base_folder = marking_win.UserData.base_folder;
[~,base_name] = fileparts(base_folder);
marking_file = [base_folder filesep base_name '_fepsp_markings.mat'];
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
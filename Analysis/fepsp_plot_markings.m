function [Visu_Win] = fepsp_plot_markings(varargin)
% Show all markings and slope areas (if avaliable) on all traces.
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "get_protocol.m" for
%                 more info.
%   markings    - Struct. Hold all of the responses "starts" & "peaks"
%                 loacations. See "fepsp_markings.m" for more info.

% INPUT (optional):
%   results     - Struct. Hold all of the analysis results & info. 
%                 See "fepsp_analyse.m" for more info. If empty slope
%                 measurement area won't be plotted.
%                 Default: [].
%   intens      - Numeric vector of stimulus intensities used for each
%                 column of traces. Typically provided in units of uA. if
%                 empty will be set as a unit increasing vector.
%                 Default: [].
%   traces_xlim - Numeric vector of 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty will use protocol default values. See
%                 "get_protocol.m" for more info. 
%                 Default: [].
%   traces_ylim - Numeric vector with 2 elements. voltage limits for
%                 displaying the trace [mV]. These y-limits remain constant
%                 throughout all intensities. If empty will be set
%                 widest values in each channel (excluding stimulus
%                 artifact). 
%                 Default: [].
%   dt          - Non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "get_protocol.m" for more info.
%                 Default: 2.
%
% OUTPUT:
%   Visu_Win    - figure array, channel (row) x intensity (column). Plots
%                 of each traces group with its markings.
%
% CALL:
%   Visu_Win = fepsp_plot_markings("traces", <cell array>, "fs", <numeric
%   scalar>, "protocol_id", <string scalar>,"markings",<struct
%   scalar>,"results",<struct scalar or empty>,"intens",<numeric
%   vector>,"traces_Xlimit",<numeric 2 elements or empty>,
%   "traces_Ylimit",<numeric 2 elements or empty>,"dt",<numeric
%   non-negative scalar>)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('traces',        [], @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',            [], @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',   [], @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('markings',      [], @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('results',       [], @(x) (isstruct(x) && isscalar(x)) || isempty(x))
p.addParameter('intens',        [], @(x) validateattributes(x,{'numeric'},{'vector'}))
p.addParameter('traces_xlim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('dt',            2,  @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))

p.parse(varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
markings        = p.Results.markings;
results         = p.Results.results;
intens          = p.Results.intens;
traces_xlim     = sort(p.Results.traces_xlim);
traces_ylim     = sort(p.Results.traces_ylim);
dt              = p.Results.dt;

% validate intens
if isempty(intens)
    intens = -(size(traces, 2) : -1 : 1);
elseif numel(intens) ~= size(traces, 2)
    error('number of intensities does not fit size of traces')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get protocol info
protocol_info = get_protocol("protocol_id",protocol_id,"fs",fs,"dt",dt);
% if user gave traces_xlim, take them over protocol default
if ~isempty(traces_xlim)
    protocol_info.traces_xlim = traces_xlim;
end
% forcing Tstamps to be row vectors solve cat error due to different
% behaviour when using a matrix or a vector as linear indexes to a vector.
protocol_info.Tstamps = protocol_info.Tstamps(:)'; 

% choose y limits
if isempty(traces_ylim)
    % if user didn't defined y limits, take the 110% of the widest needed
    % to present the all of the traces in the channel, without any removed traces
    for iChan = size(traces,1):-1:1
        all_channel.traces = [traces{iChan,:}];
        all_channel.starts = [markings.starts{iChan,:}];
        all_channel.peaks = [markings.peaks{iChan,:}];
        deleted_traces = all(isnan(all_channel.peaks),1) | all(isnan(all_channel.starts),1);
        all_channel.traces(:,deleted_traces) = nan;
        traces_ylim(iChan,:) = [min(all_channel.traces(protocol_info.response.win,:),[],'all','omitnan'),...
            max(all_channel.traces(protocol_info.response.win,:),[],'all','omitnan')]*1.1;
    end
else
    % user defined y limits - use them for all channels
    traces_ylim = repmat(traces_ylim,size(traces,1),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTraces_group = numel(traces):-1:1
    % extract data relevant to this traces group:
    loop_starts = markings.starts{iTraces_group};
    loop_peaks = markings.peaks{iTraces_group};
    loop_traces_group = traces{iTraces_group};

    % fix deleted traces to NaN. Identify by columns that are NaN in markings
    deleted_traces = all(isnan(loop_starts),1) | all(isnan(loop_peaks),1);
    loop_traces_group(:,deleted_traces) = nan;

    % convert the marks to 1 so they can be used as index. Because
    % loop_traces_group for traces that were deleted are now NaN, any value
    % plotted using them won't show.
    loop_starts(:,deleted_traces) = 1;
    loop_peaks(:,deleted_traces) = 1;
    
    % extract data of the average trace of the same traces group
    [iChan,iInten] = ind2sub(size(traces),iTraces_group);
    loop_avg_trace = mean(loop_traces_group,2,'omitnan');
    loop_starts_avg = markings.starts_avg{iTraces_group};
    loop_peaks_avg = markings.peaks_avg{iTraces_group};
    
    % open new figure
    Visu_Win(iChan,iInten) = figure();
    
    % get number of traces in this trace group
    nTrace = size(loop_traces_group,2);
    
    % plot the traces
    traces_plot = plot(protocol_info.Tstamps,loop_traces_group);
    % Show trace number on datatips
    for iTrace = 1:nTrace
        traces_plot(iTrace).Tag = num2str(iTrace);
        r = dataTipTextRow('Trace num',ones(size(traces_plot(iTrace).XData))*iTrace);
        traces_plot(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
    end
    % plot average trace
    hold on
    traces_plot(nTrace+1) = plot(protocol_info.Tstamps,loop_avg_trace,'--b','LineWidth',2);
    traces_plot(end).DataTipTemplate.DataTipRows(end) = dataTipTextRow('Avg Trace','');
    
    % create start markers:
    % extract markers so every trace will be in a different column, rows are stimulus number 
    linear_starts = sub2ind(size(loop_traces_group),loop_starts,repmat(1:nTrace,protocol_info.nStim,1));
    starts_X = [protocol_info.Tstamps(loop_starts), protocol_info.Tstamps(loop_starts_avg)'];
    start_Y = [loop_traces_group(linear_starts), loop_avg_trace(loop_starts_avg)];
    if protocol_info.nStim == 1
        % if there is only 1 trace, markers will be on a row vector,
        % resulting in 1 line obj for all markers. Place every marker on a
        % separated column, when all the other values are NaN - that way
        % only 1 marker is plotted for each line obj
        starts_X = diag(starts_X);
        starts_X(~eye(size(starts_X))) = nan;
        start_Y = diag(start_Y);
        start_Y(~eye(size(start_Y))) = nan;
    end
    % plot the start marker
    start_scatter = plot(starts_X,start_Y,'k*');
    
    % create peak markers:
    % same as start markers, see above for comments
    linear_peaks = sub2ind(size(loop_traces_group),loop_peaks,repmat(1:size(loop_traces_group,2),protocol_info.nStim,1));
    peaks_X = [protocol_info.Tstamps(loop_peaks), protocol_info.Tstamps(loop_peaks_avg)'];
    peaks_Y = [loop_traces_group(linear_peaks), loop_avg_trace(loop_peaks_avg)];
    if protocol_info.nStim == 1
        peaks_X = diag(peaks_X);
        peaks_X(~eye(size(peaks_X))) = nan;
        peaks_Y = diag(peaks_Y);
        peaks_Y(~eye(size(peaks_Y))) = nan;
    end
    peaks_scatter = plot(peaks_X,peaks_Y,'kO');
    
    % limits, labels, title & legend
    xlim(protocol_info.traces_xlim)
    ylim(traces_ylim(iChan,:))
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    legend([start_scatter(end);peaks_scatter(end)],{'response Start','response Peak'},'Location','best')
    title(sprintf('Channel %G Intensity %G',iChan,intens(iInten)))
    
    % if results area available, plot also the slope measurement area.
    slope_area_mark = gobjects(nTrace+1,protocol_info.nStim);
    if ~isempty(results)
        % extract data relevant to this traces group & average trace
        loop_slope_win = sort(results.all_traces.slope_win{iTraces_group}); % Sorted to prevent colon operator break
        loop_slope_win(isnan(loop_slope_win)) = 1; % deleted traces are already NaN - convert loop_slope_win from NaN to 1 so it can be used to index
        loop_avg_slope_win = sort(results.avg_traces.slope_win{iTraces_group});
        
        % colour each slope area
        for iStim = protocol_info.nStim:-1:1
            for iTrace = nTrace:-1:1
                % colour slope area for each trace
                slope_area_mark(iTrace,iStim) = plot(protocol_info.Tstamps(loop_slope_win(1,iStim,iTrace):loop_slope_win(2,iStim,iTrace)),...
                    loop_traces_group(loop_slope_win(1,iStim,iTrace):loop_slope_win(2,iStim,iTrace),iTrace),'r','LineWidth',3);
                % Make marked area transparent
                slope_area_mark(iTrace,iStim).Color(4) = 0.5;
                % Make marked area not appear in legend
                slope_area_mark(iTrace,iStim).Annotation.LegendInformation.IconDisplayStyle = 'off';
                % Make hover data tips take from the line under the marked area
                slope_area_mark(iTrace,iStim).PickableParts = 'none';
            end
            % colour slope area for average trace
            slope_area_mark(nTrace+1,iStim) = plot(protocol_info.Tstamps(loop_avg_slope_win(1,iStim):loop_avg_slope_win(2,iStim)),...
                loop_avg_trace(loop_avg_slope_win(1,iStim):loop_avg_slope_win(2,iStim)),'r','LineWidth',3);
            slope_area_mark(nTrace+1,iStim).Color(4) = 0.5;
            slope_area_mark(nTrace+1,iStim).Annotation.LegendInformation.IconDisplayStyle = 'off';
            slope_area_mark(nTrace+1,iStim).PickableParts = 'none';
        end
        
        % add 1 slope area colour mark to legend, from average trace
        slope_area_mark(nTrace+1,1).Annotation.LegendInformation.IconDisplayStyle = 'on';
        slope_area_mark(nTrace+1,1).DisplayName = sprintf('Slope: %.3g%% to %.3g%%',results.slope_area*100);
    end
    
    % create checkboxes to toggle each trace visibility, along with all its markers
    for iTrace = nTrace:-1:1
        trace_checkbox(iTrace) = uicontrol('String',sprintf('Trace %d',iTrace),'Style','checkbox','Value',1,...
            'Callback',@(obj,~) vis_tog(obj,traces_plot(iTrace),start_scatter(iTrace),peaks_scatter(iTrace),slope_area_mark(iTrace,:)));
    end
    % disable checkboxes for deleted traces
    [trace_checkbox(deleted_traces).Enable] = deal('off');
    
    % make checkboxes for average trace
    trace_checkbox(nTrace+1) = uicontrol('String','Avg trace','Style','checkbox','Value',1,...
        'Callback',@(obj,~) vis_tog(obj,traces_plot(nTrace+1),start_scatter(nTrace+1),peaks_scatter(nTrace+1),slope_area_mark(nTrace+1,:)));
    
    % make button for quick showing or hiding all traces
    check_all = uicontrol('String','Hide/Show all','Style','pushbutton','Tooltip','Hide / Show all','Callback',@(~,~) hide_show_all(trace_checkbox));
    
    % align traces checkbox vertically
    align(trace_checkbox,'Left','Fixed',5)
    % align button with lowest checkbox
    align([trace_checkbox(1),check_all],'Fixed',5,'Bottom')
    
    % fix controls too short for their strings
    trace_checkbox(nTrace+1).Position(3) = 70;
    check_all.Position(3) = 70;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vis_tog(obj,traces_plot,start_scatter,peak_scatter,slope_area_mark)
% hide / show visibility of a trace and every marker on it when checkbox is
% unchecked / checked.
%
% INPUT
%   obj         - uicontrol object. The checkbox that unchecked / checked,
%                 caller of this callback.
%   traces_plot - line object. The plot of a single trace connected to this
%                 checkbox.
%   start_scatter
%               - line object. Mark of response start for this trace.
%   peak_scatter
%               - line object. Mark of response peak for this trace.
%   slope_area_mark
%               - GraphicsPlaceholder or line object. If line, the colour
%                 marking of the slope measurement area for this trace.


if obj.Value == 0
    % checkbox is unchecked, hide trace & markers
    traces_plot.Visible = 'off';
    start_scatter.Visible = 'off';
    peak_scatter.Visible = 'off';
    if ~isa(slope_area_mark,'matlab.graphics.GraphicsPlaceholder')
        [slope_area_mark.Visible] = deal('off');
    end
else
    % checkbox is checked, show trace & markers
    traces_plot.Visible = 'on';
    start_scatter.Visible = 'on';
    peak_scatter.Visible = 'on';
    if ~isa(slope_area_mark,'matlab.graphics.GraphicsPlaceholder')
        [slope_area_mark.Visible] = deal('on');
    end
end
end


function hide_show_all(trace_checkbox)
% if at least 1 trace is showing, hide them all. Else, show all.
%
% INPUT
%   trace_checkbox
%               - uicontrol array. All the checkboxes for all the traces.

% check if any trace is still showing (checkbox is checked)
any_checked = any([trace_checkbox.Value]);

if any_checked
    % uncheck all checkboxes
    [trace_checkbox.Value] = deal(0);
else
    % check all checkboxes
    [trace_checkbox.Value] = deal(1);
end

% update traces to show / not by the status of their trace_checkbox
for iCheckbox = 1:numel(trace_checkbox)
    trace_checkbox(iCheckbox).Callback(trace_checkbox(iCheckbox))
end
end

% EOF
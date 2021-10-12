function [analysed_fepsp] = fepsp_show(varargin)
% Create a figure for each channel, summarising the fepsp experiment.
% Will create input-output graph if protocol has only 1 stimulus, and
% facilitation graph if it has more.
%
% INPUTS:
%   traces      - Required. 2d cell array of channel (row) x intensity
%                 (column) with each cell containing a 2d numeric mat of
%                 sample (row) x repetition (trace, column). This is
%                 similar to "Package Standard Cell Organization" (see also
%                 README), but inside each cell, rows are samples instead
%                 of stimulus number.
%                 This is the organized raw data.
%   fs          - Required. Numeric scalar. Sampling frequency [Hz].
%   protocol_id - Required. String scalar or char vector. ID of stimulation
%                 protocol used, for example "io","stp" or "custom" (for
%                 manually writing a new protocol). 
%                 See "GetProtocol" for more info.
%   Avg_marked_peaks -
%                 Required. 2d cell array of channel (row) x intensity
%                 (column), with each cell containing a numeric column
%                 vector, rows are stimulus number. This is similar to
%                 "Package Standard Cell Organization" (see also README),
%                 but inside each cell, there is only 1 column - for the
%                 average trace - instead of a column for every repetition
%                 (trace). Each value is the sample number of the response
%                 peak, for the average trace of the matching traces group
%                 in "traces".
%                 See fepsp_markings.
%   Avg_marked_starts -
%                 Required. Organized the same way as "Avg_marked_peaks"
%                 (see above). Each value is the sample number of the
%                 response starts, for the average trace of the matching
%                 traces group in "traces".
%                 See fepsp_markings
%   Results     - Required. Struct scalar. Holds all of the results & analysis info.
%                 Has many fields (see fepsp_analyse), but for this
%                 function only necessary:
%       slope_area - Numeric 2 elements between 0 & 1. Part of
%                    amplitude to measure the slope between. For example,
%                    if slope_area is [0.2 0.9] then slope will be measured
%                    in the response between the points that match 20% &
%                    90% of the amplitude.
%                    See fepsp_analyse.
%       all_traces - Struct scalar. Analysis results for all the traces in
%                    every trace group. In it the required fields:
%           Amp         - "Package Standard Cell Organization" (see
%                         README). Each value is the amplitude of matching
%                         response.
%           Slope       - "Package Standard Cell Organization" (see
%                         README). Each value is the slope of matching
%                         response.
%       avg_traces - Struct scalar. Analysis results for the average traces
%                    matching each traces group. In it the required field:
%           slope_window_avg - 
%                         2d cell array of channel (row) x intensity
%                         (column), with each cell containing a numeric
%                         mat, slope_area element (row, always size 2) x
%                         stimulus number (column). This is similar to
%                         "Package Standard Cell Organization" (see also
%                         README), but inside each cell, the rows are
%                         slope_area elements, and the columns are stimulus
%                         number (instead of being rows).
%                         The sample numbers that define the slope was
%                         measured with.
%                         See fepsp_analyse for more info.
%   intens      - Optional. Numeric vector describing the stimulus
%                 intensity used for each column of traces.
%                 Default: An increasing vector, from minus number of
%                 columns in traces up to minus 1. For example if traces
%                 has 3 columns, the default intens will be the vector:
%                 [-3,-2,-1].
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
%                 are constant - they remine the default values between
%                 channels.
%                 If empty ( [] ), take for each channel the limits needed
%                 to present the average traces, excluding stimulus
%                 artifact.
%                 See GetWindowsFromProtocol for more info.
%                 Default: Empty ( [] ).
%   dt          - Optional. Non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See GetWindowsFromProtocol for more info.
%                 Default: 2.
%   base_folder - Optional. The folder in which the whole project is
%                 managed. The name of this folder will be present in every
%                 file of this project, and if any required input is
%                 missing, we will try to load it from file named in the
%                 base_folder.
%                 Inputs "traces", "fs" & "protocol_id" will be loaded from
%                 file named "<base_folder_name>_fepsp_traces.mat".
%                 For example if base folder path is "C:\Users\Slutsky", the
%                 file has to be named "Slutsky_fepsp_traces.mat").
%                 Inputs of marks ("Avg_marked_peaks" &
%                 "Avg_marked_starts") will be loaded from file named
%                 "<base_folder_name>_fepsp_markings", for example
%                 "Slutsky_fepsp_markings".
%                 Analysis results ("Results") will be loaded from file named
%                 "<base_folder_name>_fepsp_results", for example
%                 "Slutsky_fepsp_results".
%                 In order to load variables from the file, they must have
%                 the same name as expected in inputs (for example,
%                 "traces").
%                 See README for more information.
%                 Default: MATLAB's current folder
%
% OUTPUTS:
%   analysed_fepsp -
%                 figure handles array. Handles to the summary figures, one
%                 for each channel.
%
% EXAMPLES
%   1) * Full Function Call *
%   analysed_fepsp = fepsp_show("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>,...
%   "Avg_marked_peaks",<cell array>,"Avg_marked_starts",<cell array>,...
%   "Results",<struct scalar>,...
%   "base_folder",<folder path>,"intens",<numeric vector>,"traces_Xlimit",<numeric 2 elements, or empty>,"traces_Ylimit",<numeric 2 elements, or empty>,...
%   "dt",<numeric non-negative scalar>)
%   % Filled example:
%   % This example shows function call, while explicitly defining any
%   % optional input default value. Notice that required inputs -
%   % "traces", "fs", "protocol_id", "Avg_marked_peaks",
%   % "Avg_marked_starts" & "Results" - do not have default value, 
%   % and an input example is shown here instead.
%   analysed_fepsp = fepsp_show("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom",...
%    "Avg_marked_peaks",{500,500;500,500},"Avg_marked_starts",{450,450;450,450},...
%    "Results",struct("slope_area",[0.2 0.9],...
%       "all_traces",struct("Amp",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)},"Slope",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)}),...
%       "avg_traces",struct("slope_window_avg",{[470;480],[470;480];[470;480],[470;480]})),...
%    "base_folder",pwd,"intens",[-2,-1],"traces_Xlimit",[],"traces_Ylimit",[],"dt",2)
%
%   2) * Defining only required inputs, using optional defaults by omitting *
%   % Any optional input can be set to default value by simply omitting it.
%   % You can omit all or only some of the optional inputs.
%   analysed_fepsp = fepsp_show("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>,...
%   "Avg_marked_peaks",<cell array>,"Avg_marked_starts",<cell array>,...
%   "Results",<struct scalar>)
%   % Filled example:
%   % The filled function call in example 1 is the same as simply:
%   analysed_fepsp = fepsp_show("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom",...
%    "Avg_marked_peaks",{500,500;500,500},"Avg_marked_starts",{450,450;450,450},...
%    "Results",struct("slope_area",[0.2 0.9],...
%       "all_traces",struct("Amp",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)},"Slope",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)}),...
%       "avg_traces",struct("slope_window_avg",{[470;480],[470;480];[470;480],[470;480]})))
%
%   3) * Loading missing required inputs from files in "base_folder" *
%   % If any required inputs ("traces", "fs", "protocol_id",
%   % "Avg_marked_peaks", "Avg_marked_starts", "Results") are missing, and
%   % base_folder contain a file containing them (as specified in
%   % "base_folder" input description), they will be loaded. For example
%   % the following code will produce the same results as the filled example 2:
%   traces = {ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)};
%   fs = 100;
%   protocol_id = "custom";
%   save("C:\User\Slutsky\Slutsky_fepsp_traces.mat","traces","fs","protocol_id")
%   Avg_marked_peaks = {500,500;500,500};
%   Avg_marked_starts = {450,450;450,450};
%   save("C:\User\Slutsky\Slutsky_fepsp_markings.mat","Avg_marked_peaks","Avg_marked_starts")
%   "Results",struct("slope_area",[0.2 0.9],...
%       "all_traces",struct("Amp",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)},"Slope",{ones(4,1),ones(5,1);ones(4,1),ones(5,1)}),...
%       "avg_traces",struct("slope_window_avg",{[470;480],[470;480];[470;480],[470;480]}));
%   save("C:\User\Slutsky\Slutsky_fepsp_results.mat","Results")
%   analysed_fepsp = fepsp_show("base_folder","C:\User\Slutsky")
%   % #Note1: You can use optional inputs by specifying them in this case, as always.
%   % #Note2: calling fepsp_markings without any inputs (fepsp_show() )
%   % will use the current MATLAB folder as is the default for "base_folder".

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
p.addParameter('Avg_marked_starts',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('Avg_marked_peaks',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('Results',[],@(x) validateattributes(x,{'struct'},{'scalar'}))

% Function management arguments
p.addParameter('traces_Xlimit',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('traces_Ylimit',[],@(x) isnumeric(x) && (numel(x)==2 || isempty(x)))
p.addParameter('dt',2,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))
p.addParameter('base_folder', pwd,@isfolder)


parse(p, varargin{:})

% Data arguments
traces = p.Results.traces;
fs = p.Results.fs;
protocol_id = p.Results.protocol_id;
intens = p.Results.intens;
Avg_marked_starts = p.Results.Avg_marked_starts;
Avg_marked_peaks = p.Results.Avg_marked_peaks;
Results = p.Results.Results;

% Function management arguments
base_folder = p.Results.base_folder;
traces_Xlimit = sort(p.Results.traces_Xlimit);
traces_Ylimit = sort(p.Results.traces_Ylimit);
dt = p.Results.dt;

% Check if any required input is missing
[~,base_name] = fileparts(base_folder);

% Traces:
req_inputs_traces = {'traces','fs','protocol_id'};
req_logIND_traces = ismember(req_inputs_traces,p.UsingDefaults);
if any(req_logIND_traces)
    % If there are missing inputs, try to load them from '_fepsp_traces'
    traces_file = [base_folder,filesep,base_name,'_fepsp_traces.mat'];
    if exist(traces_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs_traces{req_logIND_traces}),traces_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(traces_file,req_inputs_traces{req_logIND_traces})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(traces,{'cell'},{'2d'},'fepsp_show','traces')
        validateattributes(fs,{'numeric'},{'scalar'},'fepsp_show','fs')
        validateattributes(protocol_id,{'string','char'},{'scalartext'},'fepsp_show','protocol_id')
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''traces'',''fs'',''protocol_id''} are required.'],sprintf('%s, ',req_inputs_traces{req_logIND_traces}),traces_file)
    end
end

% Markings:
req_inputs_marks = {'Avg_marked_starts','Avg_marked_peaks'};
req_logIND_marks = ismember(req_inputs_marks,p.UsingDefaults);
if any(req_logIND_marks)
    % If there are missing inputs, try to load them from '_fepsp_markings'
    marks_file = [base_folder,filesep,base_name,'_fepsp_markings.mat'];
    if exist(marks_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs_marks{req_logIND_marks}),marks_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(marks_file,req_inputs_marks{req_logIND_marks})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(Avg_marked_starts,{'cell'},{'2d'},'fepsp_show','Avg_marked_starts')
        validateattributes(Avg_marked_peaks,{'cell'},{'2d'},'fepsp_show','Avg_marked_peaks')
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''Avg_marked_starts'',''Avg_marked_peaks''} are required.'],...
            sprintf('%s, ',req_inputs_marks{req_logIND_marks}),marks_file)
    end
end

% Analysis results:
if ismember('Results',p.UsingDefaults)
    results_file = [base_folder,filesep,base_name,'_fepsp_results.mat'];
    if exist(results_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {''Results''} from found file "%s"\n',results_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(results_file,'Results')
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        validateattributes(Results,{'struct'},{'scalar'},'fepsp_show','Results')
    else
        error(['Missing required input/s {''Results''}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''Results''} are required.'],...
            marks_file)
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

% % Windows mark important area in the trace - baseline (before first
% stimulus), display (xlimits) and response area (expected response area,
% for excluding stimulus artifacts).
% See GetWindowsFromProtocol for more info
windows_info = GetWindowsFromProtocol(protocol_info,fs,dt,traces_Xlimit);

% Get how many intensities & stimulations there are
nIntens = length(intens);
nStim = protocol_info.nstim;

% Find the order of intensities
[intens_sorted,intens_order] = sort(intens,'ascend');

% Create time stamps
Tstamps = -protocol_info.StimLatency:(1/fs)*1000:(protocol_info.TraceRecordingLenght-protocol_info.StimLatency);
Tstamps = Tstamps';

% Get parameters about slope measured area
slope_window_avg = Results.avg_traces.slope_window_avg;
slope_area_label = sprintf('%.3g%% to %.3g%%',Results.slope_area*100);

% Get calculated slope & amplitude
Amp = Results.all_traces.Amp;
Slope =  Results.all_traces.Slope;

%% Create fepsp summary graphs

for iChan = size(traces,1):-1:1
    
    % Open figure & subplots
    analysed_fepsp(iChan) = figure();
    sgtitle(sprintf('%s - Channel #%d', base_name, iChan),'Interpreter','none')
    if nStim > 1
        subplot(length(intens)+1, 1, 1)
    else
        subplot(2, 2, [1 2])
    end
    
    % Get average traces
    for iInten = nIntens:-1:1
        % Convert deleted traces to NaN before averaging, by finding what
        % Amp columns are NaN
        deleted_traces = all(isnan(Amp{iChan,iInten}),1);
        loop_traces_group = traces{iChan,iInten};
        loop_traces_group(:,deleted_traces) = nan;
        % Create average trace
        avg_traces(:,iInten) = mean(loop_traces_group,2,'omitnan');
    end
    
    % Plot average traces and add intensity to datatips
    avg_traces_plot = plot(Tstamps,avg_traces,'LineWidth',1);
    for iInten = 1:nIntens
        avg_traces_plot(iInten).Tag = num2str(intens(iInten));
        r = dataTipTextRow('intensity',ones(size(avg_traces_plot(iInten).XData))*intens(iInten));
        avg_traces_plot(iInten).DataTipTemplate.DataTipRows(end+1) = r;
    end
    hold on
    
    % Plot starts markers:
    % Get linear_starts
    loop_starts_avg = [Avg_marked_starts{iChan,:}];
    linear_starts = sub2ind(size(avg_traces),loop_starts_avg,repmat(1:nIntens,nStim,1));
    % Plot starts via scatter
    starts_scatter = plot(Tstamps(loop_starts_avg),avg_traces(linear_starts),'k*');
    
    % Plot peaks markers:
    % Get linear_peaks
    loop_peaks_avg = [Avg_marked_peaks{iChan,:}];
    linear_ends = sub2ind(size(avg_traces),loop_peaks_avg,repmat(1:nIntens,nStim,1));
    % Plot peaks via scatter
    ends_scatter = plot(Tstamps(loop_peaks_avg),avg_traces(linear_ends),'kO');
    
    % Color slope area in each trace
    for iStim = nStim:-1:1
        for iInten = nIntens:-1:1
            % Get slope window & plot on top of samples matching it
            loop_slope_window = sort(slope_window_avg{iChan,iInten}(:,iStim)); % Sorted to prevent colon operator break
            slope_area_mark(iInten) = plot(Tstamps(loop_slope_window(1):loop_slope_window(2)),...
                avg_traces(loop_slope_window(1):loop_slope_window(2),iInten),'b','LineWidth',3);
            % Make marked area transparent
            slope_area_mark(iInten).Color(4) = 0.5;
            % Make marked area not appear in legend
            slope_area_mark(iInten).Annotation.LegendInformation.IconDisplayStyle = 'off';
            % Make hover data tips take from the line under the marked area
            slope_area_mark(iInten).PickableParts = 'none';
        end
    end
    
    % finishing touch:
    % Fix axes limit for visibility by windows
    xlim(Tstamps(windows_info.DispWindowBase))
    if isempty(traces_Ylimit)
        ylim([min(avg_traces(windows_info.AnalyseWin,:),[],'all') max(avg_traces(windows_info.AnalyseWin,:),[],'all')].*1.1)
    else
        ylim(traces_Ylimit)
    end
    % Add legend
    legend([avg_traces_plot(intens_order);starts_scatter(1);ends_scatter(1)],[string(intens_sorted),{'Measure Start','Measure End'}],...
        'Location','southeast','NumColumns',2)
    % Make 1 of the marked area appear in legend, to represent all of them
    slope_area_mark(1).Annotation.LegendInformation.IconDisplayStyle = 'on';
    slope_area_mark(1).DisplayName = ['Slope: ' slope_area_label];
    % axis labels & title
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    title('Average Traces')
    
    if nStim > 1
        % multi stimuli protocol - show facilitation.
        % For each intensity create 1 box plot for amplitude & 1 for
        % slope
        for iInten = 1:nIntens
            loop_subplot = subplot(nIntens+1, 1, iInten+1);
            
            % Normalised amplitude & slope
            loop_amp_norm = Amp{iChan, iInten}./Amp{iChan, iInten}(1,:);
            loop_slope_norm = Slope{iChan, iInten}./Slope{iChan, iInten}(1,:);
            
            % Order data so amplitude & slope are next to each other in
            % every stimulus number
            cat_data = nan(size(traces{iChan, iInten},2),2*nStim);
            col2fill = mat2cell(1:2*nStim,1,ones(1,nStim)*2);
            for iStim = nStim:-1:1
                cat_data(:,col2fill{iStim}) = [loop_amp_norm(iStim,:)' loop_slope_norm(iStim,:)'];
            end
            
            % Plot boxplots
            boxplot(cat_data(:,3:end));
            
            % Color boxplots by type
            colors = repmat({'b','r'},1,nStim);
            box_edges = findobj(loop_subplot,'Tag','Box');
            for iBox = (nStim*2-2):-1:1
                color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),colors{iBox},'FaceAlpha',.5,'PickableParts','none');
            end
            
            % finishing touch:
            % Fix ticks to be centred between boxes - make it seem like
            % the number refer both of the boxes above it
            xticks(loop_subplot,1.5:2:(nStim*2-2));
            xticklabels(string(2:nStim))
            % axis labels & title       
            xlabel('Number of Stimulation [#]')
            ylabel({'Mean Part of' 'First Stimulation' '[stim/(stim num 1)]'})
            legend(color_box(1:2),['Slope: ' slope_area_label],'Amplitude','Location','best')
        end
    else
        % Only 1 stimulation in protocol - show input-output
        
        % Get how many traces are in each intens
        nTraces = cellfun(@(x) size(x,2),traces(iChan,:));
        intens_group = repelem(string(intens),nTraces);
        
        % plot input-output for amplitude
        loop_subplot = subplot(2, 2, 3);
        Amp_all_intens = [Amp{iChan,:}];
        boxplot(loop_subplot,Amp_all_intens,intens_group,'GroupOrder',string(intens_sorted))
        xticklabels(string(intens_sorted))
        xlabel('Intensity [uA]')
        ylabel('Amplidute [mV]')
        title('Input/Output (Amplidute)')
        
        % Color boxplots 2 match avg_traces_plot colors
        box_edges = findobj(loop_subplot,'Tag','Box');
        for iBox = nIntens:-1:1
            iInten = intens_order(end-iBox+1); % Intensity is inversed & sorted between box & input intens order
            color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),avg_traces_plot(iInten).Color,'FaceAlpha',.5,'PickableParts','none');
        end
        
        % plot input-output for slope
        loop_subplot = subplot(2, 2, 4);
        slope_all_intens = [Slope{iChan,:}];
        boxplot(loop_subplot,slope_all_intens,intens_group,'GroupOrder',string(intens_sorted))
        xlabel('Intensity [uA]')
        ylabel('Slope [mV/ms]')
        title(['Input/Output (Slope: ' slope_area_label ')'])
        
        % Color boxplots 2 match avg_traces_plot colors
        box_edges = findobj(loop_subplot,'Tag','Box');
        for iBox = nIntens:-1:1
            iInten = intens_order(end-iBox+1); % Intensity is inversed & sorted between box & input intens order
            color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),avg_traces_plot(iInten).Color,'FaceAlpha',.5,'PickableParts','none');
        end
    end
end

end

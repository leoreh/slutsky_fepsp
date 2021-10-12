function [Results]= fepsp_analyse(varargin)
% Calculate amplitude and slope for each trace (repetition), using marking
% of the start & peak of each response.
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
%   marked_peaks -
%                 Required. "Package Standard Cell Organization" (see README).
%                 Each value is the sample number of the response peak.
%                 Note: Any traces removed during marking, their marked
%                 peaks will be NaN for all stimulus numbers (rows).
%                 They will be skipped during analysis.
%   marked_starts - 
%                 Required. "Package Standard Cell Organization" (see README).
%                 Each value is the sample number of the response start.
%                 Note: Any traces removed during marking, their marked
%                 starts will be NaN for all stimulus numbers (rows).
%                 They will be skipped during analysis.
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
%   Avg_marked_starts -
%                 Required. Organized the same way as "Avg_marked_peaks"
%                 (see above). Each value is the sample number of the
%                 response starts, for the average trace of the matching
%                 traces group in "traces".
%   base_folder - Optional. The folder in which the whole project is
%                 managed. The name of this folder will be present in every
%                 file of this project, and if any required input is
%                 missing, we will try to load it from file named in the
%                 base_folder.
%                 Inputs "traces", "fs" & "protocol_id" will be loaded from
%                 file named "<base_folder_name>_fepsp_traces.mat".
%                 For example if base folder path is "C:\Users\Slutsky", the
%                 file has to be named "Slutsky_fepsp_traces.mat").
%                 Inputs of marks ("marked_peaks", "marked_starts",
%                 "Avg_marked_peaks" & "Avg_marked_starts") will be loaded
%                 from file named "<base_folder_name>_fepsp_markings", for
%                 example "Slutsky_fepsp_markings".
%                 In order to load variables from the file, they must have
%                 the same name as expected in inputs (for example,
%                 "traces").
%                 See README for more information.
%                 Default: MATLAB's current folder
%   save_var    - Optional. Logical flag. Saving results at base_folder as
%                 mat file, named "<base_folder_name>_fepsp_results". For
%                 example if base path is "C:\Users\Slutsky", the saved
%                 file will be "Slutsky__fepsp_results". Will overwrite
%                 if needed.
%                 Default: true
%   slope_area  - Optional. Numeric 2 elements between 0 & 1. Part of
%                 amplitude to measure the slope between. For example, if
%                 slope_area is [0.2 0.9] then slope will be measured in
%                 the response between the points that match 20% & 90% of
%                 the amplitude. Order does not matter - slope between 20%
%                 & 90% of amplitude is the same as slope between 90% & 20%.
%                 Default: [0.2 0.9]
% OUTPUTS:
%   Results     - Struct scalar. Holds all of the results & analysis info.
%                 Has fields:
%       slope_area - The same as input slope_area, see above. Between what
%                    parts of response amplitude slope was measured.
%       all_traces - Struct scalar. Analysis results for all the traces in
%                    every trace group. Has fields:
%           Amp         - "Package Standard Cell Organization" (see
%                         README). Each value is the amplitude of matching
%                         response.
%           Slope       - "Package Standard Cell Organization" (see
%                         README). Each value is the slope of matching
%                         response.
%           slope_err   - 2d cell array of channel (row) x intensity
%                         (column), with each cell containing a 2d struct
%                         array of stimulus number (row) x repetition (trace,
%                         column). This is similar to "Package Standard
%                         Cell Organization" (see README), but a struct
%                         array inside each cell instead of numeric.
%                         Error estimation structure generated by polyfit,
%                         for each calculated slope.
%                         See polyfit for more info.
%           slope_window - 
%                         2d cell array of channel (row) x intensity
%                         (column), with each cell containing a 3d numeric
%                         mat, slope_area element (dim 1, always size 2) x
%                         stimulus number (dim 2) x repetition (trace, dim 3).
%                         This is similar to "Package Standard Cell
%                         Organization" (see README), but inside each cell
%                         there is an added dimension - the matching
%                         slope_area element, creating a 3d numeric array.
%                         The sample number that define the window. Each
%                         element matches the corresponding element in
%                         slope_area. For example if slope are is [0.5
%                         0.1], the first element in slope_window is the
%                         sample that match 50% of amplitude inside of
%                         response, and the second is the sample number of
%                         10%.
%       avg_traces - Struct scalar. Analysis results for the average traces
%                    matching each traces group. Has the same fields as
%                    all_traces, but ending with "_avg":
%           Amp_avg     - 2d cell array of channel (row) x intensity
%                         (column), with each cell containing a numeric
%                         column vector, rows are stimulus number. This is
%                         similar to "Package Standard Cell Organization"
%                         (see also README), but inside each cell, there is
%                         only 1 column - for the average trace - instead
%                         of a column for every repetition (trace).
%                         Each value is the amplitude of matching
%                         response, in the average trace.
%           slope_avg   - Orginazed the same way as Amp_avg (see above). 
%                         Each value is the slope of matching response, in
%                         the average trace. 
%           slope_err   - 2d cell array of channel (row) x intensity
%                         (column), with each cell containing a column vector struct
%                         rows are stimulus number. This is similar to
%                         "Package Standard Cell Organization" (see
%                         README), but a struct array inside each cell
%                         instead of numeric, and only 1 column (for the
%                         average trace) instead of a column for each trace.
%                         Error estimation structure generated by polyfit,
%                         for each calculated slope.
%                         See polyfit for more info.
%           slope_window_avg - 
%                         2d cell array of channel (row) x intensity
%                         (column), with each cell containing a numeric
%                         mat, slope_area element (row, always size 2) x
%                         stimulus number (column). This is similar to
%                         "Package Standard Cell Organization" (see also
%                         README), but inside each cell, the rows are
%                         slope_area elements, and the columns are stimulus
%                         number (instead of being rows).
%                         Same as field "slope_window" in "all_traces", but
%                         for the average trace.
% EXAMPLES
%   1) * Full Function Call *
%   Results = fepsp_analyse("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>,...
%       "marked_peaks",<cell array>,"marked_starts",<cell array>,"Avg_marked_peaks",<cell array>,"Avg_marked_starts",<cell array>,...
%       "base_folder",<folder path>,"save_var",<logical flag>,"slope_area",<numeric 2 elements, 0<= & >=1>)
%   % Filled example:
%   % This example shows function call, while explicitly defining any
%   % optional input default value. Notice that required inputs -
%   % "traces", "fs", "protocol_id", "marked_peaks", "marked_starts",
%   % "Avg_marked_peaks" & "Avg_marked_starts" - do not have default value,
%   % and an input example is shown here instead.
%   Results = fepsp_analyse("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom",...
%       "marked_peaks",{ones(1,4)*500,ones(1,5)*500;ones(1,4)*500,ones(1,4)*500},"marked_starts",{ones(1,4)*450,ones(1,5)*450;ones(1,4)*450,ones(1,4)*450},...
%       "Avg_marked_peaks",{500,500;500,500},"Avg_marked_starts",{450,450;450,450},...
%       "base_folder",pwd,"save_var",true,"slope_area",[0.2 0.9])
%
%   2) * Defining only required inputs, using optional defaults by omitting *
%   % Any optional input can be set to default value by simply omitting it.
%   % You can omit all or only some of the optional inputs.
%   Results = fepsp_analyse("traces",<cell array>,"fs",<numeric scalar>,"protocol_id",<string scalar>,...
%       "marked_peaks",<cell array>,"marked_starts",<cell array>,"Avg_marked_peaks",<cell array>,"Avg_marked_starts",<cell array>)
%   % Filled example:
%   % The filled function call in example 1 is the same as simply:
%   Results = fepsp_analyse("traces",{ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)},"fs",100,"protocol_id","custom",...
%       "marked_peaks",{ones(1,4)*500,ones(1,5)*500;ones(1,4)*500,ones(1,4)*500},"marked_starts",{ones(1,4)*450,ones(1,5)*450;ones(1,4)*450,ones(1,4)*450},...
%       "Avg_marked_peaks",{500,500;500,500},"Avg_marked_starts",{450,450;450,450})
%
%   3) * Loading missing required inputs from files in "base_folder" *
%   % If any required inputs ("traces", "fs", "protocol_id",
%   % "marked_peaks", "marked_starts", "Avg_marked_peaks" &
%   % "Avg_marked_starts") are missing, and base_folder contain a file
%   % containing them (as specified in "base_folder" input description),
%   % they will be loaded.
%   % For example the following code will produce the same results as
%   % the filled example 2:
%   traces = {ones(1000,4),ones(1000,5);ones(1000,4),ones(1000,5)};
%   fs = 100;
%   protocol_id = "custom";
%   save("C:\User\Slutsky\Slutsky_fepsp_traces.mat","traces","fs","protocol_id")
%   marked_peaks = {ones(1,4)*500,ones(1,5)*500;ones(1,4)*500,ones(1,4)*500};
%   marked_starts = {ones(1,4)*450,ones(1,5)*450;ones(1,4)*450,ones(1,4)*450};
%   Avg_marked_peaks = {500,500;500,500};
%   Avg_marked_starts = {450,450;450,450};
%   save("C:\User\Slutsky\Slutsky_fepsp_markings.mat","marked_peaks","marked_starts","Avg_marked_peaks","Avg_marked_starts")
%   Results = fepsp_analyse("base_folder","C:\User\Slutsky")
%   % #Note1: You can use optional inputs by specifying them in this case, as always.
%   % #Note2: calling fepsp_analyse without any inputs (fepsp_markings() )
%   % will use the current MATLAB folder as is the default for "base_folder".
%
%   4) * Changing slope area *
%   % Slope area is defined between start & peak marks, between 2 points
%   % that match certien part of the total response amplitude. Taking only
%   % part of the response "wave" this way reduce the effect that
%   % non-linear wave parts have on the slope. You can change what area you
%   % are using for slope mesuarment by using the slope_area input:
%   Results = fepsp_analyse(___,"slope_area",<numeric 2 elements, 0<= & >=1>)
%   % You can do that with any other syntax.
%   % Filled example: (simillar to example 3, everything is presaved)
%   Results = fepsp_analyse("base_folder","C:\User\Slutsky","slope_area",[0.2 0.5])


% Note1: Change polyfit to fit. Problem - too few points case.
% Note2: Path conflict with the fepsp_analysis in slutskycode.

%% Input Parser - Define & validate All Inputs

%Give Option to enter inputs via struct, instead of name-value
p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

% Traces
p.addParameter('traces',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',[],@(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',[],@(x) validateattributes(x,{'string','char'},{'scalartext'}))

% Markings
p.addParameter('marked_starts',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('marked_peaks',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('Avg_marked_starts',[],@(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('Avg_marked_peaks',[],@(x) validateattributes(x,{'cell'},{'2d'}))

% Function management arguments
p.addParameter('base_folder',pwd,@isfolder)
p.addParameter('save_var',true,@(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.addParameter('slope_area',[0.20 0.90],@(x) validateattributes(x,'numeric',{'numel',2,'>=',0,'<=',1}))

p.parse(varargin{:})

% Traces
traces = p.Results.traces;
fs = p.Results.fs;
protocol_id = p.Results.protocol_id;

% Markings
marked_peaks = p.Results.marked_peaks;
marked_starts = p.Results.marked_starts;
Avg_marked_starts = p.Results.Avg_marked_starts;
Avg_marked_peaks = p.Results.Avg_marked_peaks;

% Function management arguments
base_folder = p.Results.base_folder;
save_var = p.Results.save_var;
slope_area = p.Results.slope_area(:)'; % Make sure it is row vec

% Check if any required input is missing
% Traces:
req_inputs_traces = {'traces','fs','protocol_id'};
req_logIND_traces = ismember(req_inputs_traces,p.UsingDefaults);
if any(req_logIND_traces)
    % If there are missing inputs, try to load them from '_fepsp_traces'
    [~,base_name] = fileparts(base_folder);
    traces_file = [base_folder filesep base_name '_fepsp_traces.mat'];
    if exist(traces_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs_traces{req_logIND_traces}),traces_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(traces_file,req_inputs_traces{req_logIND_traces})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(traces,{'cell'},{'2d'},'fepsp_analyse','traces')
        validateattributes(fs,{'numeric'},{'scalar'},'fepsp_analyse','fs')
        validateattributes(protocol_id,{'string','char'},{'scalartext'},'fepsp_analyse','protocol_id')
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''traces'',''fs'',''protocol_id''} are required.'],sprintf('%s, ',req_inputs_traces{req_logIND_traces}),traces_file)
    end
end

% Markings:
req_inputs_marks = {'marked_peaks','marked_starts','Avg_marked_starts','Avg_marked_peaks'};
req_logIND_marks = ismember(req_inputs_marks,p.UsingDefaults);
if any(req_logIND_marks)
    % If there are missing inputs, try to load them from '_fepsp_markings'
    [~,base_name] = fileparts(base_folder);
    marks_file = [base_folder filesep base_name '_fepsp_markings.mat'];
    if exist(marks_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs_marks{req_logIND_marks}),marks_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(marks_file,req_inputs_marks{req_logIND_marks})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(marked_peaks,{'cell'},{'2d'},'fepsp_analyse','marked_peaks')
        validateattributes(marked_starts,{'cell'},{'2d'},'fepsp_analyse','marked_starts')
        validateattributes(Avg_marked_starts,{'cell'},{'2d'},'fepsp_analyse','Avg_marked_starts')
        validateattributes(Avg_marked_peaks,{'cell'},{'2d'},'fepsp_analyse','Avg_marked_peaks')
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''marked_peaks'',''marked_starts'',''Avg_marked_starts'',''Avg_marked_peaks''} are required.'],...
            sprintf('%s, ',req_inputs_marks{req_logIND_marks}),marks_file)
    end
end


%% Crete helper variables
% Get protocol struct info from protocol_id
protocol_info = GetProtocol(protocol_id);

% Tstamps matching each trace
Tstamps = -protocol_info.StimLatency:(1/fs)*1000:(protocol_info.TraceRecordingLenght-protocol_info.StimLatency);
Tstamps = Tstamps';


%% Initiate Analysed Outputs

% Non Average Traces Outputs:
% cell electrode-by-intensity - inside numeric stim-by-trace:
[Amp,Slope] = deal(cellfun(@(x) nan(protocol_info.nstim,size(x,2)),traces,'UniformOutput',0));
% cell electrode-by-intensity - inside struct array stim-by-trace:
fit_err_struct = struct('R',[],'df',[],'normr',[]);
slope_err = cellfun(@(x) repmat(fit_err_struct,protocol_info.nstim,size(x,2)),traces,'UniformOutput',0);
% cell electrode-by-intensity - inside numeric 3d: point-by-stim-by-trace
% - sample first & last sample of the window used to calculate slope, as
% defined by the percent of amplitude in "slope_area"
slope_window = cellfun(@(x) nan(2,protocol_info.nstim,size(x,2)),traces,'UniformOutput',0);

% Average Traces Outputs:
% cell electrode-by-intensity - inside numeric stim (rows) col vec:
[Amp_avg,slope_avg] = deal(repmat({nan(protocol_info.nstim,1)},size(traces)));
% cell electrode-by-intensity - inside struct array stim (rows) col vec:
slope_err_avg = repmat({repmat(fit_err_struct,protocol_info.nstim,1)},size(traces));
% cell electrode-by-intensity - inside numeric 2d: point-by-stim 
% - sample first & last sample of the window used to calculate slope, as
% defined by the percent of amplitude in "slope_area"
slope_window_avg = repmat({nan(2,protocol_info.nstim)},size(traces));


%% Analyse Stage

% Polynomial not unique might occur when using polyfit. As we are using
% polyfit multiple times, it will be uninformative and clutter the command
% window. We will silence it, and present a clear warning instead when it
% rises (see later)
warning('off','MATLAB:polyfit:PolyNotUnique')

for iTraces_group = numel(traces):-1:1 %Run through each traces group
    %% Calculate Amplitude for all the traces
    
    % Extract data relevant to this traces group:
    % Of all the traces of this trace group
    loop_starts = marked_starts{iTraces_group};
    loop_peaks = marked_peaks{iTraces_group};
    loop_traces_group = traces{iTraces_group};
    % Fix deleted traces to NaN. Identify by columns that are NaN in markings
    deleted_traces = all(isnan(loop_starts),1) | all(isnan(loop_peaks),1);
    loop_traces_group(:,deleted_traces) = nan;
    % Convert the marks to 1 so they can be used as index. Because
    % loop_traces_group for traces that were deleted are now NaN, amplitude
    % will also be NaN.
    loop_starts(:,deleted_traces) = 1;
    loop_peaks(:,deleted_traces) = 1;
    
    % Of the average trace of the same traces group
    [iChan,iInten] = ind2sub(size(traces),iTraces_group); % the subscript of this traces group match the channel & intensity
    loop_avg_trace = mean(loop_traces_group,2,'omitnan'); % Create avg traces from loop_traces_group to take into accuont deletion
    loop_starts_avg = Avg_marked_starts{iTraces_group};
    loop_peaks_avg = Avg_marked_peaks{iTraces_group};
    
    % Calculate amplitude for all the traces:
    % First, convert indices from subscripts to linear, in order to
    % maintain 1 index refer 1 sample.
    linear_starts = sub2ind(size(loop_traces_group),loop_starts,repmat(1:size(loop_traces_group,2),protocol_info.nstim,1));
    linear_ends = sub2ind(size(loop_traces_group),loop_peaks,repmat(1:size(loop_traces_group,2),protocol_info.nstim,1));
    % Amplitude is simply absolute of: start of response minus end of response
    Amp{iTraces_group} = abs(loop_traces_group(linear_starts) - loop_traces_group(linear_ends));
    
    % Redo for average trace of the same traces group
    Amp_avg{iTraces_group} = abs(loop_avg_trace(loop_starts_avg) - loop_avg_trace(loop_peaks_avg));
    %% Calculate each trace slopes
    % We will loop through every stimulation (and its matching response, as a result) given inside each protocol repetition (trace). In
    % each trace and in the average trace, we will calculate the slope. See
    % subfunction "calculate_slope" for slope calculation algorithm.
    
    for iStim = protocol_info.nstim:-1:1
        
        % Find slopes for each trace in traces group
        for iTrace = size(loop_traces_group,2):-1:1
            % Skip deleted traces - their amplitude will be NaN
            if isnan(Amp{iTraces_group}(iStim,iTrace))
                fprintf('Skiped deleted trace %d in channel %d & intensity %d\n',iTrace,iChan,iInten)
                continue
            end
            
            % Create a window in which the response occur, between
            % response start & peak
            responce_win = loop_starts(iStim,iTrace):loop_peaks(iStim,iTrace);
            
            % Calculate Slope:
            [slope_window{iTraces_group}(:,iStim,iTrace),Slope{iTraces_group}(iStim,iTrace),slope_err{iTraces_group}(iStim,iTrace)] = ...
                calculate_slope(Tstamps,loop_traces_group(:,iTrace),responce_win,Amp{iTraces_group}(iStim,iTrace),slope_area);
            
            % There are a few cases when there are not enough points
            % between percentiles - usually when response amplitude is very
            % small. We silenced those warnings earlier, and we will
            % produce one informative warning instead
            if strcmp(lastwarn,'Polynomial is not unique; degree >= number of data points.')
                lastwarn('') % Clear lastwarn to prevent warning reappear next loop
                warning(['Its seems there are too few points when searching slope, in Channel %G, Intens %G, StimNum %G in trace %G. '...
                    'As a results, relatively big fit error may occur. '...
                    'This may be due to Start & Peak marker too close, combination of very small slope & big "max_time_tol", or a biphasic event between Start & Peak markers. ' ...
                    'If Slope is expected to be close to 0 in this trace, you can ignore. '...
                    'Check slope_err_10_50 & slope_err_20_90, at the {%G,%G}(%G,%G) For more info about the fit error'],...
                    iChan,iInten,iStim,iTrace,...
                    iChan,iInten,iStim,iTrace)
            end
        end
        
        % Redo for average trace of the same traces group:
        
        % Create a window in which the response occur, from user
        % define end & start.
        responce_win = loop_starts_avg(iStim):loop_peaks_avg(iStim);
        
        % Calculate Slope:
        [slope_window_avg{iTraces_group}(:,iStim),slope_avg{iTraces_group}(iStim),slope_err_avg{iTraces_group}(iStim)] = ...
            calculate_slope(Tstamps,loop_avg_trace,responce_win,Amp_avg{iTraces_group}(iStim),slope_area);
        
        % There are a few cases when there are not enough points
        % between percentiles - usually when response amplitude is very
        % small. We silenced those warnings earlier, and we will produce
        % one informative warning instead.
        if strcmp(lastwarn,'Polynomial is not unique; degree >= number of data points.')
            lastwarn('')
            warning(['Its seems there are too few points when searching slope, in Channel %G, Intens %G, StimNum %G in the average trace. '...
                'As a results, relatively big fit error may occur. '...
                'This may be due to Start & Peak marker too close, combination of very small slope & big "max_time_tol", or a biphasic event between Start & Peak markers. ' ...
                'If Slope is expected to be close to 0 in this trace, you can ignore. '...
                'Check slope_err_10_50_avg & slope_err_20_90_avg, at the {%G,%G}(%G) For more info about the fit error'],...
                iChan,iInten,iStim,...
                iChan,iInten,iStim)
        end
    end
end
% Now we finished using polyfit - lets restore the warning to its regular
% "on" status, to avoid missing it when it is needed.
warning('on','MATLAB:polyfit:PolyNotUnique')

%% Organize outputs & Save

% Pack all outputs in struct
Results.all_traces = struct('Amp',{Amp},'Slope', {Slope},'slope_err',{slope_err},'slope_window',{slope_window});
Results.avg_traces = struct('Amp_avg',{Amp_avg},'slope_avg',{slope_avg},'slope_err_avg',{slope_err_avg},'slope_window_avg',{slope_window_avg});
Results.slope_area = slope_area;
if save_var
    % Create file name from base_folder
    [~,base_name] = fileparts(base_folder);
    results_file = [base_folder filesep base_name '_fepsp_results.mat'];
    fprintf('Saving results as "%s"\n',results_file)
    save(results_file,'Results')
end
end


function [slope_window,Slope,slope_err] = calculate_slope(Tstamps,trace,responce_win,Amp,slope_area)
% Calculate slope in a single trace between 2 point inside of responce_win
% that match slope_area part of amplitude.
% Algorithm: Calculate all of the amplitudes of the points in
% responce_win (absolute difference from first point in responce_win). Find
% the points whose amplitude value is the closest to the part of Amp
% specified by slope_area. Fit a line to all the points between the two
% points found - this line slope is the slope.
% INPUTS
%   Tstamps - Numeric vector. time [ms] that match the trace.
%   trace   - Numeric vector. single data sample repetition in [mV].
%   responce_win -
%             Positive ascending integer vector. Samples number for the
%             response (start to peak), to be taken from trace & Tstamps
%   Amp     - Numeric scalar, calculated response amplitude for this
%             trace [mV].
%   slope_area - 
%             Numeric 2 elements between 0 & 1. Part of amplitude to
%             measure the slope between. For example, if slope_area is [0.2
%             0.9] then slope will be measured in the response between the
%             points that match 20% & 90% of the amplitude. Order does not
%             matter - slope between 20% & 90% of amplitude is the same as
%             slope between 90% & 20%.
% OUTPUT
%   slope_window - 
%             Numeric 2 elements positive integer.
%             The sample number that define the window. Each element matches
%             the corresponding element in slope_area. For example if slope
%             are is [0.5 0.1], the first element in slope_window is the
%             sample that match 50% of amplitude inside of response, and the
%             second is the sample number of 10%.
%   Slope   - Numeric scalar, calculated response slope for this trace
%             [mV/ms].
%   slope_err - 
%             Struct scalar. Error estimation structure generated by
%             polyfit. See polyfit for more info.

% First find all the response’s points amplitude - absolute
% difference from response start data sample
points_amp = abs(trace(responce_win(1)) - trace(responce_win));
% Create a vector of "how much are slope_area's percentiles of calculated
% amplitude?"
p_vec = Amp*slope_area;
% Using implicit expansion, find data samples with the closest
% amplitude to wanted percentile of general amplitude
[~,slope_window] = min(abs(points_amp-p_vec),[],1);
% Right now, slope_window match the responce_win instead of sample numbers
% in the general trace. Convert them.
slope_window = slope_window + responce_win(1) - 1;
% If there are both a significant biphasic event inside the
% response, the lower percentage might be found after the
% higher one. Will break colon operator later. So sort them.
slope_window_sorted = sort(slope_window,'ascend');
% Notice that we keep output unsorted, so the first the 1st element in
% slope_window matching the first percentile in slope_area

% Calculate Slope by fitting a line. Save results & fit error to output:
[fit_params,slope_err] = ...
    polyfit(Tstamps(slope_window_sorted(1):slope_window_sorted(2)), trace(slope_window_sorted(1):slope_window_sorted(2)), 1);
% Polyfit first output is in descending power order. We don't
% need the intercept (power == 0) value.
Slope = fit_params(1);
end


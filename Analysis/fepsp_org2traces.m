function [traces] = fepsp_org2traces(varargin)
% Organize data from whole channel into traces.
% Takes data samples (data_samples), Sampling Frequency (fs), Stimulation Protocol
% (protocol_id) & location of each protocol repetition first stimuli
% (stim_locs). Returns traces (protocol repetitions) Organized so that they
% can be passed to "fepsp_markings".
% Algorithm: Assume that data samples are continues & that all channels were
% recorded at the same time. For each stimulation in stim_locs, take all
% the points matching the entire stimulation protocol, and extract them
% from all channels. Finally, organize all the extracted traces in similar
% way to a 'Package Standard Cell Organization' (that is explained in the
% README).
%
% INPUT
%   data_samples - 
%                 Required. 2d numeric mat of sample (row) x channel
%                 (column). Data recorded during the experiment. 
%                 This function assume that the data was recorded
%                 continuously (according to the sampling frequency) - sample
%                 2 was taken one sampling period after sample 1 and so on.
%   fs          - Required. Numeric scalar. Sampling frequency [Hz].
%   protocol_id - Required. String scalar or char vector. ID of stimulation
%                 protocol used, for example "io","stp" or "custom" (for
%                 manually writing in a new protocol).
%                 See "GetProtocol" for more info.
%   stim_locs   - Required. Cell vector, in each cell there is a numeric
%                 row vector of positive integers, without any NaN or inf.
%                 Each cell is a different intensity. Each cell contains
%                 the sample numbers matching the first stimulus given in
%                 each protocol repetition (traces).
%                 Those stimulation are assumed to be recorded in all
%                 the channels simultaneously.
%   base_folder - Optional. Path to the folder in which the whole project is
%                 managed. The name of this folder will be present in every
%                 file of this project, and if any required input is
%                 missing, we will try to load it from files in the
%                 base_folder. In order to load from file, it must be named
%                 "<base_folder_name>_fepsp_base_data.mat"
%                 (for example if base folder is named "Slutsky", the
%                 file must be named "Slutsky_fepsp_base_data.mat").
%                 The file must contain the required inputs arguments:
%                 "data_samples","fs","protocol_id","stim_locs".
%                 See README for more information.
%                 Default: MATLAB's current folder.
%   remove_DC   - Optional. Logical flag. Remove DC component.
%                 Default: true
%   save_var    - Optional. Logical flag. Saving results at base_folder as
%                 mat file, named "<base_folder_name>_fepsp_traces". For
%                 example if base path is "C:\Users\Slutsky", the saved
%                 file will be "Slutsky_fepsp_traces". Will overwrite
%                 if needed.
%                 Default: true
%
% OUTPUT
%   traces      - 2d cell array of channel (row) x intensity (column) with
%                 each cell containing a 2d numeric mat of sample (row) x
%                 repetition (trace, column). This is similar to "Package Standard
%                 Cell Organization" (see also README), but inside each
%                 cell, rows are samples instead of stimulus
%                 number.
%                 Those are "traces" - protocol repetitions, packed in
%                 "traces groups" ordered by channel & intensity. Other
%                 package functions expect this organization.
%
% EXAMPLES
%   1) * Full Function Call *
%   traces = ...
%   fepsp_org2traces("data_samples",<numeric mat>,"fs",<numeric scalar>,"protocol_id",<string scalar>,"stim_locs",<cell vec>,...
%   "base_folder",<folder path>,"remove_DC",<logical flag>,"save_var",<logical flag>)
%   % Filled example:
%   % This example shows function call, while explicitly defining any
%   % optional input default value. Notice that required inputs -
%   % "data_samples","fs", "protocol_id" & "stim_locs" - do not have
%   % default value, and an input example is shown here instead.
%   traces = ...
%   fepsp_org2traces("data_samples",[1 5;2 6;3 7;4 8],"fs",2,"protocol_id","custom","stim_locs",{[2],[3 4]},...
%   "base_folder",pwd,"remove_DC",false,"save_var",false)
%
%   2) * Defining only required inputs, using optional defaults by omitting *
%   % Any optional input can be set to default value by simply omitting it.
%   % You can omit all or only some of the optional inputs.
%   traces = ...
%   fepsp_org2traces("data_samples",<numeric mat>,"fs",<numeric scalar>,"protocol_id",<string scalar>,"stim_locs",<cell vec>)
%   % Filled example:
%   % The filled function call in example 1 is the same as simply:
%   traces = ...
%   fepsp_org2traces("data_samples",[1 5;2 6;3 7;4 8],"fs",2,"protocol_id","custom","stim_locs",{[2],[3 4]})
%
%   3) * Loading missing required inputs from files in "base_folder" *
%   % If any required inputs ("traces","fs","protocol_id') are missing, and
%   % base_folder contain a file containing them (as specified in
%   % "base_folder" input description), they will be loaded. 
%   % For example the following code will produce the same results as
%   % the filled example 2:
%   data_samples = [1 5;2 6;3 7;4 8]
%   fs = 2;
%   protocol_id = "custom";
%   stim_locs = {[2],[3 4]};
%   save("C:\User\Slutsky\Slutsky_fepsp_base_data.mat","data_samples","fs","protocol_id","stim_locs")
%   traces = fepsp_org2traces("base_folder","C:\User\Slutsky")
%   % #Note1: You can use optional inputs by specifying them in this case, as always.
%   % #Note2: calling fepsp_org2traces without any inputs (fepsp_org2traces() )
%   % will use the current MATLAB folder, as is the default for "base_folder".
%
% TROUBLESHOOTING
%   Data too big for MATLAB array (Out of Memory error):
%       1) Use memmapfile + taking advantage of Copy-On-Write MATLAB
%       behaviour: (see https://www.mathworks.com/help/matlab/matlab_prog/avoid-unnecessary-copies-of-data.html)
%           % First Call memmapfile to map the file, then pass all the data
%           % to the function. For example:
%               data_file_name = 'path_to_data_samples_File.dat';
%               data_samples_map = memmapfile(data_file_name,'Format',{'double',[nSamp nCh],'data_samples');
%               traces = ...
%                   fepsp_org2traces("data_samples",data_samples_map.Data.data_samples,"fs",<numeric scalar>,"protocol_id",<string scalar>,"stim_locs",<cell vec>)
%           % Since the function does not make any changes to data_samples,
%           % it continues to read them from the file as it copies.
%       2) Cut your data to shorter segments. Match stim_locs to each
%       segment. Verify that you leave enough space in each segment after
%       the last stim_locs & before the first stim_locs for the whole trace
%       to fit.
%
%   I have the time the stimulation was given (relative to start of the
%   recording), not the sample number:
%       Convert them to sample number. This can be easily done by
%       multiplying them by multiplying the time in seconds by sampling
%       frequency. If the first sample in data_samples correspond to time
%       0, add 1 after multiplying. Round the result to compensate for any
%       equipment inaccuracies.
%       For example, if stimulation was given
%       30 millisec after recording start, first sample match time 0,
%       and sampling frequency is 20000, then it matches sample number 601:
%       ((30/1000)*20000) + 1 = 601.
% 
%       Code example for transforming:
%       % Convert stimulation locations from times to sample numbers. Assume
%       % that "stim_times" is a cell vector similar to "stim_locs", but each
%       % stimulation location is in time from start of recording (which started
%       % at time 0) instead of data samples:
%       stim_locs = cellfun(@(x) round(x*fs+1),stim_times);
%       fepsp_org2traces("data_samples",<numeric mat>,"fs",<numeric scalar>,"protocol_id",<string scalar>,"stim_locs",stim_locs)
%
%   See also fepsp_load, fepsp_markings

% Note1 - removeDC need to include slope removal - find slope between the
%         first & last points of the baseline & remove the created line
%         from the whole trace. If signal is 0.5X+15, remove the line 0.5X
%         (X is time or sample number).
% Note2 - Should this function also include dt? Can we add dt to protocol?
%         if we can't, dt passes all functions.

%% Input Parser - Define & validate All Inputs
p = inputParser;

% Give Option to enter inputs via struct, instead of name-value
p.StructExpand = true;
p.KeepUnmatched = true;

% Data arguments
p.addParameter('data_samples',[],@(x) validateattributes(x,{'numeric'},{'2d'}))
p.addParameter('fs',[],@(x) validateattributes(x,{'numeric'},{'scalar','positive'}))
p.addParameter('protocol_id','',@(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('stim_locs',[],@(x) validateattributes(x,{'cell'},{'vector'}))

% Function management arguments
p.addParameter('base_folder',pwd,@isfolder)
p.addParameter('remove_DC',true,@(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.addParameter('save_var',true,@(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.parse(varargin{:})

% Data arguments
data_samples   =   p.Results.data_samples;
fs             =   p.Results.fs;
protocol_id    =   p.Results.protocol_id;
stim_locs      =   p.Results.stim_locs;

% Function management arguments
remove_DC      =   logical(p.Results.remove_DC);
base_folder    =   p.Results.base_folder;
save_var       =   p.Results.save_var;

% Check if any required input is missing
req_inputs = {'data_samples','fs','protocol_id','stim_locs'};
req_logIND = ismember(req_inputs,p.UsingDefaults);
if any(req_logIND)
    % If there are missing inputs, try to load them from '_fepsp_traces'
    [~,base_name] = fileparts(base_folder);
    base_data_file = [base_folder filesep base_name '_fepsp_base_data.mat'];
    if exist(base_data_file, 'file')
        % File exist. Notify user about the loading, and load by varnames.
        % If any variable needed to be loaded is missing, error instead of warning.
        fprintf('Loading missing required input/s {%s} from found file "%s"\n',sprintf('%s, ',req_inputs{req_logIND}),base_data_file)
        
        warning('error','MATLAB:load:variableNotFound')
        load(base_data_file,req_inputs{req_logIND})
        warning('on','MATLAB:load:variableNotFound')
        
        % Validate loaded variables, the same way as when given as inputs.
        % Any variable given as input will pass this check, as it has passed already.
        validateattributes(data_samples,{'numeric'},{'2d'},'fepsp_org2traces','data_samples')
        validateattributes(fs,{'numeric'},{'scalar'},'fepsp_org2traces','fs')
        validateattributes(protocol_id,{'string','char'},{'scalartext'},'fepsp_org2traces','protocol_id')
        validateattributes(stim_locs,{'cell'},{'vector'},'fepsp_org2traces','stim_locs')
        
    else
        error(['Missing required input/s {%s}. Was unable to find file "%s" to load them from.',...
            ' Inputs {''data_samples'',''fs'',''protocol_id'',''stim_locs''} are required.'],sprintf('%s, ',req_inputs{req_logIND}),base_data_file)
    end
end

%% Create Helpers Values

% Get how many intensities, channels & samples there are
nCh = size(data_samples,2);
nSamp = size(data_samples,1);
nIntens = size(stim_locs,2);

% Get protocol struct info from protocol_id
protocol_info = GetProtocol(protocol_id);

% Base for cutting latter - IND matching the IND of each trace that
% will be cut, centred at the first stimulus in the protocol.
% Floor because if there wasn't enough time to take another sample, it wasn't taken
first_sample = -floor(protocol_info.StimLatency.*fs/1000);
last_sample = floor(protocol_info.TraceRecordingLenght.*fs/1000) + first_sample;
trace_mold = first_sample:1:last_sample;

%% Cut all channels by the stimulation times given

% Create IND matrix for the indices matching each trace,
% centred on the first stimulus in the protocol. The matrix is
% created using Implicit expansion (trace_mold transposed is
% column vector, ALL_stim_locs is row vector) and as a result in it each
% column is a trace (different stimulus), and each row is index of
% different data sample.
stim_locs_cat = [stim_locs{:}];
traces_unorg = trace_mold' + stim_locs_cat;


try
    % Extract trace from all the channels. MATLAB subscript order is by
    % column, which match traces_base shape.
    all_traces = data_samples(traces_unorg,:);
catch err
    % In case any of the traces have out of bound indices -
    % meaning protocol recording should end after channel recording
    % stopped, or protocol should start before recording started -
    % bring a nice error message. Else, unexpected error occurred -
    % just rethrow it
    Overflow = max(traces_unorg,[],1) > nSamp;
    Underflow = min(traces_unorg,[],1) < 1;
    if ~any(Overflow) && ~any(Underflow)
        rethrow(err)
    end
    
    % Find in which intensities cells are matching which problem-causing stim_locs
    nStim = cellfun(@numel,stim_locs);
    nStim_culm = cumsum(nStim);
    [~,intens_overflow] = max(nStim_culm'./find(Overflow) >= 1,[],1);
    [~,intens_underflow] = max(nStim_culm'./find(Underflow) >= 1,[],1);
    
    % Generate nice error message, describing which stim_loc are causing
    % which type of error:
    err_msg = '';
    if any(Overflow)
        over_msg = sprintf(' Stimulus locations {%s} in corresponding intensities cells {%s} produce traces that end after recording finished.',...
            sprintf('%d,',stim_locs_cat(Overflow)),sprintf('%d,',intens_overflow));
        err_msg = [err_msg over_msg];
    end
    if any(Underflow)
        under_msg = sprintf(' Stimulus locations {%s} in corresponding intensities cells {%s} produce traces that start before recording started.',...
            sprintf('%d,',stim_locs_cat(Underflow)),sprintf('%d,',intens_underflow));
        err_msg = [err_msg under_msg];
    end
    error([err_msg ' Check if stim_locs & Protocol are correct, and if data has the expected number of data samples.'])
end

%% Organize traces in a Package Standard Cell Organization

% Separate stimulation repetitions by intensity, by "reverting" 
% the stim_locs horzcat done earlier
nStim = cellfun(@numel,stim_locs);
traces = mat2cell(all_traces,nStim*length(trace_mold),ones(nCh,1))';

% Separate in each intensity between the different repetitions (traces).
for iIntens = 1:nIntens
    for iChannel = 1:nCh
        traces{iChannel,iIntens} = reshape(traces{iChannel,iIntens},length(trace_mold),[]);
    end
end

% Remove DC - calculate DC from pre-stimulus period
if remove_DC
    baseline = [1 find(trace_mold == -floor(fs/1000))]; % From signal start, to 1 ms before stimulation (in order to ignore artifact)
    for iIntens = 1:nIntens
        for iChannel = 1:nCh
            traces{iChannel,iIntens} = rmDC(traces{iChannel,iIntens}, 'dim', 1, 'win', baseline);
        end
    end
end

if save_var
    % Create file name from base_folder
    [~,base_name] = fileparts(base_folder);
    traces_file = [base_folder filesep base_name '_fepsp_traces.mat'];
    fprintf('Saving traces as "%s"\n',traces_file)
    save(traces_file,'traces','fs','protocol_id')
end
end
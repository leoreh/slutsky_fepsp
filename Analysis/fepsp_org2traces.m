function [traces] = fepsp_org2traces(varargin)

% organizes continuous data (data_in) into a cell array of traces for
% downstream analysis. Requires indices of stimulation onset (stim_locs)
% arranged in a cell array according to intensities. Organization is done
% according to a predefined stimulation protocol. 
% 
% For more information see https://github.com/leoreh/slutsky_fepsp.
%
% INPUT (required):
%   data_in     - 2d numeric mat of sample (row) x channel
%                 (column). Data recorded during the experiment. 
%   fs          - numeric scalar. Sampling frequency of the
%                 recorded data [Hz]. 
%   protocol_id - string or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "get_protocol.m" for
%                 more info.
%   stim_locs   - 1d cell array of intensities with each cell
%                 containing a numeric row vector of stimulus onset indices
%                 as samples. e.g. stim_locs{3} = [5000 1000] means that
%                 two stimulations at the 3rd intensity were given at
%                 sample 5000 and 10000.

% INPUT (optional):
%   base_path   - string or char. Full path to where the output should be
%                 saved. The name of the last folder in base_path will be
%                 the prefix for all saved data. e.g. base_path =
%                 'lh85_211012_132000' than the output of fepsp_org2traces
%                 will be 'lh85_211012_132000_fepsp_traces'.
%                 Default: pwd.
%   rmv_trend   - logical flag. Remove median from data_in.
%                 Default: true
%   save_var    - logical flag. Save output in base_path as
%                 mat file named "<base_name>_fepsp_traces".
%                 Default: true
%
% OUTPUT:
%   traces      - 2d cell array of channel (row) x intensity (column) with
%                 each cell containing a 2d numeric mat of sample (row) x
%                 repetition (column). This is similar to "standard cell"
%                 (see README) but inside each cell, rows are samples
%                 instead of stimulus number.
%
% CALL:
%   traces = fepsp_org2traces("data_in", <numeric mat>, "fs", <numeric
%   scalar>, "protocol_id", <string scalar>, "stim_locs", <cell vec>,
%   "base_path", <folder path>, "rmv_trend", <logical flag>, "save_var",
%   <logical flag>);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('data_in',       [],     @(x) validateattributes(x, {'numeric'}, {'2d'}))
p.addParameter('fs',            [],     @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
p.addParameter('protocol_id',   '',     @(x) validateattributes(x, {'string', 'char'}, {'scalartext'}))
p.addParameter('stim_locs',     [],     @(x) validateattributes(x, {'cell'}, {'vector'}))
p.addParameter('base_path',     pwd,    @isfolder)
p.addParameter('rmv_trend',     true,   @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))
p.addParameter('save_var',      true,   @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))
p.parse(varargin{:})

data_in         =   p.Results.data_in;
fs              =   p.Results.fs;
protocol_id     =   p.Results.protocol_id;
stim_locs       =   p.Results.stim_locs;
rmv_trend       =   logical(p.Results.rmv_trend);
base_path       =   p.Results.base_path;
save_var        =   p.Results.save_var;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data params
nCh = size(data_in, 2);
nSamp = size(data_in, 1);
nIntens = size(stim_locs, 2);
nStim = cellfun(@numel, stim_locs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice recording according to stim indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create single trace template in samples according to protocol
protocol_info = get_protocol("protocol_id",protocol_id,"fs",fs);
trace_temp = floor(protocol_info.Tstamps*fs/1000)';

% snip traces from data by creating a mat of indices 
traces_matIdx = trace_temp' + [stim_locs{:}];

try
    traces_mat = data_in(traces_matIdx,:);

catch err
    % check if error is due to traces out of bound
    overflow = max(traces_matIdx, [], 1) > nSamp;
    underflow = min(traces_matIdx, [], 1) < 1;
    if ~any(overflow) && ~any(underflow)
        rethrow(err)
    end
    
    % find problematic intensities
    nStim_culm = cumsum(nStim);
    [~, intens_overflow] = max(nStim_culm' ./ find(overflow) >= 1, [], 1);
    [~, intens_underflow] = max(nStim_culm' ./ find(underflow) >= 1, [], 1);
    
    err_msg = '';
    if any(overflow)
        over_msg = sprintf(' trace {%s} in intensity {%s} ends after recording.',...
            sprintf('%d,', stim_locs_cat(overflow)), sprintf('%d,', intens_overflow));
        err_msg = [err_msg over_msg];
    end
    if any(underflow)
        under_msg = sprintf(' trace {%s} in intensity {%s} starts before recording.',...
            sprintf('%d,', stim_locs_cat(underflow)), sprintf('%d,', intens_underflow));
        err_msg = [err_msg under_msg];
    end
    error([err_msg ' Check stim_locs, protocol and length of data.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize in cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create cell array
traces = mat2cell(traces_mat, nStim * length(trace_temp), ones(nCh, 1))';

% reshape each cell 
for iIntens = 1 : nIntens
    for iCh = 1 : nCh
        traces{iCh, iIntens} = reshape(traces{iCh, iIntens}, length(trace_temp), []);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detrend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rmv_trend
    % build regressor to find the trend of a trace using the first and last 5
    % ms. This is to exclude the evoked response and stimulus artifact from the
    % linear fit. Based on matrix form of polynomial regression.
    lip = 5 * fs / 1000;
    reg_len = 2 * lip;
    trace_len = length(trace_temp);
    
    use_idx = [1 : lip, trace_len : -1 : trace_len - lip + 1];
    w = [0 : reg_len - 1] / (reg_len - 1);
    w = [reshape(w, reg_len, []), ones(reg_len, 1)];
    [q, r] = qr(w, 0);
    w_full = [0 : trace_len - 1] / (trace_len - 1);
    w_full = [reshape(w_full, trace_len, []), ones(trace_len, 1)];
    
    for iIntens = 1 : nIntens
        for iCh = 1 : nCh
            trend_coeff = (r \ q' * traces{iCh, iIntens}(use_idx, :));
            trend{iCh, iIntens} = w_full * trend_coeff;
            x{iCh, iIntens} = traces{iCh, iIntens} - trend{iCh, iIntens};
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if save_var
    [~, base_name] = fileparts(base_path);
    traces_file = fullfile(base_path, [base_name '_fepsp_traces.mat']);
    save(traces_file, 'traces', 'fs', 'protocol_id')
end

end

% EOF
function [analysed_fepsp] = fepsp_show(varargin)
% Create a figure for each channel, summarising the fepsp experiment.
% Will create input-output graph if protocol has only 1 stimulus, and
% facilitation graph if it has more.
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "get_protocol.m" for
%                 more info.
%   markings    - Struct. Hold all of the responses "starts" & "peaks"
%                 loacations. See "fepsp_markings.m" for more info.
%   results     - Struct. Hold all of the analysis results & info. 
%                 See "fepsp_analyse.m" for more info.

% INPUT (optional):
%   intens      - numeric vector of stimulus intensities used for each
%                 column of traces. Typically provided in units of uA. if
%                 empty will be set as a unit increasing vector.
%                 Default: [].
%   traces_xlim - numeric vector of 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty will use protocol default values. See
%                 "get_protocol.m" for more info. 
%                 Default: [].
%   traces_ylim - numeric vector with 2 elements. voltage limits for
%                 displaying the trace [mV]. These y-limits remain constant
%                 throughout all intensities. If empty will be set
%                 widest values in each channel (excluding stimulus
%                 artifact). 
%                 Default: [].
%   dt          - non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "get_protocol.m" for more info.
%                 Default: 2.
%
% OUTPUT:
%   analysed_fepsp
%               - figure handles array. Handles to the summary figures, one
%                 for each channel.
%
% CALL:
%   analysed_fepsp = fepsp_show("traces", <cell array>, "fs", <numeric
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
p.addParameter('results',       [], @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('intens',        [], @(x) (isnumeric(x) && isvector(x)) || isempty(x))
p.addParameter('traces_xlim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('dt',            2,  @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))

parse(p, varargin{:})

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
    intens = -(size(traces,2):-1:1);
elseif numel(intens) ~= size(traces,2)
    error('Number of inten does not match number of intensities in traces (number of columns in the cell array)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get protocol info
protocol_info = fepsp_getProtocol("protocol_id",protocol_id,"fs",fs,"dt",dt);
% if user gave traces_xlim, take them over protocol default
if ~isempty(traces_xlim)
    protocol_info.traces_xlim = traces_xlim;
end

% get how many intensities & stimulations there are
nIntens = length(intens);
nStim = protocol_info.nStim;

% find the order of intensities
[intens_sorted,intens_order] = sort(intens,'ascend');

% get parameters about slope measured area
slope_window_avg = results.avg_traces.slope_win;
slope_area_label = sprintf('%.3g%% to %.3g%%',results.slope_area*100);

% get calculated slope & amplitude
Amp = results.all_traces.Amp;
Slope =  results.all_traces.Slope;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iChan = size(traces,1):-1:1
 
    % open figure & subplots
    analysed_fepsp(iChan) = figure();
    sgtitle(sprintf('Channel #%d', iChan),'Interpreter','none')
    if nStim > 1
        subplot(length(intens)+1, 1, 1)
    else
        subplot(2, 2, [1 2])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot average traces for all intensities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get average traces
    for iInten = nIntens:-1:1
        % convert deleted traces to NaN before averaging, by finding what
        % Amp columns are NaN
        deleted_traces = all(isnan(Amp{iChan,iInten}),1);
        loop_traces_group = traces{iChan,iInten};
        loop_traces_group(:,deleted_traces) = nan;
        % create average trace
        avg_traces(:,iInten) = mean(loop_traces_group,2,'omitnan');
    end
    
    % plot average traces and add intensity to datatips
    avg_traces_plot = plot(protocol_info.Tstamps,avg_traces,'LineWidth',1);
    for iInten = 1:nIntens
        avg_traces_plot(iInten).Tag = num2str(intens(iInten));
        r = dataTipTextRow('intensity',ones(size(avg_traces_plot(iInten).XData))*intens(iInten));
        avg_traces_plot(iInten).DataTipTemplate.DataTipRows(end+1) = r;
    end
    hold on
    
    % create start markers:
    % convert indices from subscripts to linear, in order to maintain 1
    % index refer 1 sample.
    loop_starts_avg = [markings.starts_avg{iChan,:}];
    linear_starts = sub2ind(size(avg_traces),loop_starts_avg,repmat(1:nIntens,nStim,1));
    % plot starts as scatter
    starts_scatter = plot(protocol_info.Tstamps(loop_starts_avg),avg_traces(linear_starts),'k*');
    
    % plot peaks markers:
    % get linear index
    loop_peaks_avg = [markings.peaks_avg{iChan,:}];
    linear_peaks = sub2ind(size(avg_traces),loop_peaks_avg,repmat(1:nIntens,nStim,1));
    % plot peaks as scatter
    peaks_scatter = plot(protocol_info.Tstamps(loop_peaks_avg),avg_traces(linear_peaks),'kO');
    
    % colour slope area in each trace
    for iStim = nStim:-1:1
        for iInten = nIntens:-1:1
            % get slope window & plot on top of samples matching it
            loop_slope_window = sort(slope_window_avg{iChan,iInten}(:,iStim)); % sorted to prevent colon operator break
            slope_area_mark(iInten) = plot(protocol_info.Tstamps(loop_slope_window(1):loop_slope_window(2)),...
                avg_traces(loop_slope_window(1):loop_slope_window(2),iInten),'r','LineWidth',3);
            % make marked area transparent
            slope_area_mark(iInten).Color(4) = 0.5;
            % make marked area not appear in legend
            slope_area_mark(iInten).Annotation.LegendInformation.IconDisplayStyle = 'off';
            % make hover data tips take from the line under the marked area
            slope_area_mark(iInten).PickableParts = 'none';
        end
    end
    
    % finishing touch:
    % fix axes limit for visibility by windows
    xlim(protocol_info.traces_xlim)
    if isempty(traces_ylim)
        ylim([min(avg_traces(protocol_info.response.win,:),[],'all') max(avg_traces(protocol_info.response.win,:),[],'all')].*1.1)
    else
        ylim(traces_ylim)
    end
    
    % add legend
    legend([avg_traces_plot(intens_order);starts_scatter(1);peaks_scatter(1)],[string(intens_sorted),{'response Start','response Peak'}],...
        'Location','southeast','NumColumns',2)
    
    % make 1 of the marked area appear in legend, to represent all of them
    slope_area_mark(1).Annotation.LegendInformation.IconDisplayStyle = 'on';
    slope_area_mark(1).DisplayName = ['Slope: ' slope_area_label];
    % axis labels & title
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    title('Average Traces')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot analysis results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nStim > 1
        % multi stimuli protocol - show facilitation.
        % for each intensity create 1 box plot for amplitude & 1 for slope
        for iInten = 1:nIntens
            loop_subplot = subplot(nIntens+1, 1, iInten+1);
            
            % normalised amplitude & slope
            loop_amp_norm = Amp{iChan, iInten}./Amp{iChan, iInten}(1,:);
            loop_slope_norm = Slope{iChan, iInten}./Slope{iChan, iInten}(1,:);
            
            % order data so amplitude & slope are next to each other in
            % every stimulus number
            cat_data = nan(size(traces{iChan, iInten},2),2*nStim);
            col2fill = mat2cell(1:2*nStim,1,ones(1,nStim)*2);
            for iStim = nStim:-1:1
                cat_data(:,col2fill{iStim}) = [loop_amp_norm(iStim,:)' loop_slope_norm(iStim,:)'];
            end
            
            % plot boxplots
            boxplot(cat_data(:,3:end));
            
            % colour boxplots by type
            colors = repmat({'b','r'},1,nStim);
            box_edges = findobj(loop_subplot,'Tag','Box');
            for iBox = (nStim*2-2):-1:1
                color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),colors{iBox},'FaceAlpha',.5,'PickableParts','none');
            end
            
            % finishing touch:
            % fix ticks to be centred between boxes - make it seem like
            % the number refer both of the boxes above it
            xticks(loop_subplot,1.5:2:(nStim*2-2));
            xticklabels(string(2:nStim))
            % axis labels & title       
            xlabel('Number of Stimulation [#]')
            ylabel({'Mean Part of' 'First Stimulation' '[stim/(stim num 1)]'})
            legend(color_box(1:2),['Slope: ' slope_area_label],'Amplitude','Location','best')
        end
    else
        % only 1 stimulation in protocol - show input-output
        
        % get how many traces are in each intens
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
        
        % color boxplots 2 match avg_traces_plot colors
        box_edges = findobj(loop_subplot,'Tag','Box');
        for iBox = nIntens:-1:1
            iInten = intens_order(end-iBox+1); % intensity is inversed & sorted between box & input intens order
            color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),avg_traces_plot(iInten).Color,'FaceAlpha',.5,'PickableParts','none');
        end
        
        % plot input-output for slope
        loop_subplot = subplot(2, 2, 4);
        slope_all_intens = [Slope{iChan,:}];
        boxplot(loop_subplot,slope_all_intens,intens_group,'GroupOrder',string(intens_sorted))
        xlabel('Intensity [uA]')
        ylabel('Slope [mV/ms]')
        title(['Input/Output (Slope: ' slope_area_label ')'])
        
        % color boxplots 2 match avg_traces_plot colors
        box_edges = findobj(loop_subplot,'Tag','Box');
        for iBox = nIntens:-1:1
            iInten = intens_order(end-iBox+1); % intensity is inversed & sorted between box & input intens order
            color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),avg_traces_plot(iInten).Color,'FaceAlpha',.5,'PickableParts','none');
        end
    end
end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example data is a ~4 min recording from a mouse implanted with 3
% electrodes directed at CA1 stratum radiadum (ch1 and ch2) and startum
% pyramidale (ch 3). recording is saved in flat binary format (parameters
% listed below). if you need help loading data from any other format (e.g.
% winWCP, etc.) please contact the technical author.

% change basepath to the full path of the test data on your pc. this should
% be the only manual adjustment necassary to run this script.
basepath = 'D:\Code\slutsky_fepsp\Analysis\ExampleData';

% recording params
fs = 20000;                 % sampling frequency
nChans = 3;                 % number of channels in file
precision = 'int16';        % data saved as int16
b2mv = 0.195 / 1000;        % conversion from bits to microvolts
intens = [20 : 20 : 100];   % stim intensities used during recording [uA]
protocol_id = 'io';         % input/output stimulus protocol

% load data 
cd(basepath)
[~, basename] = fileparts(basepath);
filename = fullfile(basepath, [basename, '.dat']);
info = dir(filename);
nSamps = info.bytes / 2 / nChans;
fid = fopen(filename, 'r');
precisionTxt = sprintf('%d*%s=>double', nChans, precision);
data_in = [fread(fid, [nChans nSamps], precisionTxt)]' * b2mv;
fclose(fid);

% load sample indices of stimulus onset arranged in a cell according to the
% different intensities
load([basename, '_stimLocs.mat'])

% graphics ----------------------------------------------------------------
fh = figure;
fh = fepsp_graphics(fh);                        % set graphics
xstamps = [1 : nSamps] / fs;                    % recording duration [s]
alphaIdx = linspace(0.3, 1, length(intens));    % color tranparency 
for ich = 1 : nChans
    sb(ich) = subplot(nChans, 1, ich);

    plot(xstamps, data_in(:, ich), 'k', 'LineWidth', 0.5)
    axis tight
    hold on
    for iintens = 1 : length(intens)
        stimVec = stim_locs{iintens} / fs;
        ph = plot([stimVec; stimVec], ylim, '--', 'LineWidth', 2,...
            'Color', [1, 0, 0, alphaIdx(iintens)]);     
    end
    xlabel('Time [s]')
    ylabel('Voltage [mV]')
    title(['Channel ' num2str(ich)])
    sb(ich) = fepsp_graphics(sb(ich));          % set graphics
end
linkaxes(sb, 'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fepsp analysis pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step 1        organizes the raw data in the cell format
traces          = fepsp_org2traces('data_in', data_in, 'fs', fs,...
                'protocol_id', protocol_id, 'stim_locs', stim_locs);

% step 2        opens a gui for marking the traces
marking_win     = fepsp_markings("traces", traces, "fs", fs,...
                "protocol_id", protocol_id,"base_path", basepath,...
                "intens", intens, "traces_Xlimit", [], "traces_Ylimit", [],...
                "dt", 2, "max_jitter", 0.5, "fast_mark", true);

% step 2.5      wait until you finished marking (closed the window) & 
%               load the markings and updated traces
waitfor(marking_win)
load([basename, '_fepsp_markings.mat'], "markings")
load([basename, '_fepsp_traces.mat'], "traces")

% step 3        analyzes traces according to the manual markings
results         = fepsp_analyse("traces", traces, "fs", fs,...
                "protocol_id", protocol_id, "markings", markings,...
                "base_path", basepath, "save_var", true, "slope_area", [0.2 0.9]);

% step 4        dispalys the results
analysed_fepsp  = fepsp_summaryPlot("traces", traces, "fs", fs,...
                "protocol_id", protocol_id, "markings", markings, "results", results,...
                "base_path", basepath, "intens", intens);


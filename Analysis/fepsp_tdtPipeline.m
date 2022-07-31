function fepsp_tdtPipeline(varargin)

% organizes recordings from tdt. currently as wrapper. main problem is that
% tdt sucks! for some blocks the length of two stores (e.g. Raw1 and stim)
% defer even if they were sampled at the same frequency. this appears
% random (i.e. does not depend on block duration and can go either way. the
% workaround here is to shorten / elongate stim according to Raw1. this
% assumes the missing / additional samples are at the end of the block and
% are not important (so far has proven to be true). currently designed to
% handle the creation of a new session folder. this requires an .xml file
% in mousepath with mousename. 
% 
% INPUT
%   basepath        char. dir with .wcp files
%   blocks          numeric. see tdt2dat
%   mapch           numeric. see tdt2dat
%   rmvch           numeric. see tdt2dat
%   ch              numeric. index to channel in data stream with the fepsp
%                   signal. before rmv and map channels
%   store           char. see tdt2dat {'Raw1'}
%   protocol_id  char. can be 'freerun', 'io', or 'stp'
%   recsuffix       char. suffix to add to basename 
%   intens          numeric. intensity values of stimulations [uA]
% 
% DEPENDENCIES
%   slutsky_fepsp functions
% 
% 25 apr 22 LH
% 
% TO DO LIST


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'blocks', []);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'ch', 1, @isnumeric);
addOptional(p, 'store', 'Raw1', @ischar);
addOptional(p, 'protocol_id', 'io', @ischar);
addOptional(p, 'recsuffix', '', @ischar);
addOptional(p, 'intens', [], @isnumeric);

parse(p, varargin{:})
basepath         = p.Results.basepath;
blocks           = p.Results.blocks;
mapch            = p.Results.mapch;
rmvch            = p.Results.rmvch;
ch               = p.Results.ch;
store            = p.Results.store;
protocol_id      = p.Results.protocol_id;
recsuffix        = p.Results.recsuffix;
intens           = p.Results.intens;

if isempty(intens)
    intens = 1 : length(wcpfiles);
end

% validate input
if length(intens) ~= length(blocks)
    warning('check intens and blocks')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load from tank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)

% get fs from block header
blockpath = fullfile(basepath, ['block-', num2str(blocks(1))]);
blockHead = TDTbin2mat(blockpath, 'HEADERS', 1);
fsStim = blockHead.stores.Stim.fs;
fs = blockHead.stores.(store).fs;

% create dat
% datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
%     'chunksize', 300, 'mapch', mapch', 'rmvch', rmvch, 'clip', {});

% load dat and stim, and clip if necessary
blockfiles = dir('block*');
blocknames = {blockfiles.name};
blocknames = natsort(blocknames);
stim = []; sig = [];
for iblock = 1 : length(blocks)
    blockpath = fullfile(basepath, blocknames{blocks(iblock)});
    dat = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
        'T1', 0, 'T2', 0, 'CHANNEL', ch);
    dat = dat.streams.(store).data;
    s = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', 'Stim',...
        'T1', 0, 'T2', 0);
    s = s.streams.Stim.data;
    
    difflen(iblock) = length(s) - length(dat);
    if difflen(iblock) ~= 0
        warning('length stim and raw not equal')
        if difflen(iblock) > 0
            s(end : -1 : end - difflen(iblock) + 1) = [];
        else
            s(length(s) : length(dat)) = 0;
        end
    end
    stim = [stim, s];
    sig = [sig, dat];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize stim indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert stim binary to indices
stimIdx = find(stim);
stimIdx = [stimIdx(1), stimIdx(find(diff(stimIdx) > fsStim) + 1)];

% downsample 
fsRat = fs / fsStim;
stimIdx = round(stimIdx * fsRat);

% organize stims according to intensity
nstims = 45;
nintens = length(intens);
stimLocs = num2cell(reshape(1 : nstims, nstims / nintens, nintens), 1);
for iintens = 1 : nintens
    stimLocs{iintens} = stimIdx(stimLocs{iintens});
end

% organize in cell according to blocks
% csamps = round(cumsum(datInfo.nsamps) * fsRat);
% csamps = [1, csamps];
% stimLocs = cell(length(blocks), 1);
% for iblock = 1 : length(blocks)
%     blockIdx = stimIdx > csamps(iblock) & stimIdx < csamps(iblock + 1);
%     stimLocs{iblock} = stimIdx(blockIdx);
% end

% debugging ---------------------------------------------------------------
% fh = figure;
% plot([1 : length(sig)] / fs, sig)
% hold on
% plot([stimIdx; stimIdx] / fs, ylim, '--k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize session folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% move to folder
[~, basename] = fileparts(basepath);
% fnames = dir(['*' basename '*']);
% recname = [basename, '_', recsuffix];
% newpath = fullfile(basepath, recname);
% mkdir(newpath)
% for ifile = 1 : length(fnames)
%     if ~fnames(ifile).isdir
%         newname = strrep(fnames(ifile).name, basename, recname);
%         newfile = fullfile(newpath, newname);
%         movefile(fnames(ifile).name, newfile, 'f')
%     end
% end
newpath = basepath;
recname = basename;

% get xml from mousepath
% mousepath = fileparts(basepath);
% [~, mousename] = fileparts(mousepath);
% xmlfile = dir(fullfile(mousepath, '*xml'));
% newname = strrep(xmlfile.name, mousename, recname);
% newfile = fullfile(newpath, newname);
% copyfile(fullfile(mousepath, xmlfile.name), newfile)

% start working on session
% cd(newpath)
% basepath = pwd;
% basename = recname;
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fepsp pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMatlabGraphics(true)

% step 1        organizes the raw data in the cell format
traces          = fepsp_org2traces('data_in', double(sig') / 1000, 'fs', fs,...
                'protocol_id', protocol_id, 'stim_locs', stimLocs);

% step 2        opens a gui for marking the traces
marking_win     = fepsp_markings("traces", traces, "fs", fs,...
                "protocol_id", protocol_id,"base_path", basepath,...
                "intens", intens, "traces_Xlimit", [], "traces_Ylimit", [],...
                "dt", 2, "max_jitter", 0.5, "fast_mark", false);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add info to results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stim times in relation to entire recording
load([basename, '.', store, '.datInfo.mat'])
if blocks > 1
    addsamps = cumsum(datInfo.nsamps(1 : blocks - 1));
else
    addsamps = 0;
end
stimIdx = (stimIdx + addsamps(end)) / fs;

resultsfile = [basename, '_fepsp_results.mat'];
load(resultsfile, 'results')
results.info.tdt.ch = ch;
results.info.tdt.difflen = difflen;
results.info.intens = intens;
results.info.stimLocs = stimLocs;
results.info.stimIdx = stimIdx;
save(resultsfile, 'results')

end 

% EOF
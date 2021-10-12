
% prep
protocol_id = 'io';
intens = IntensIO;
basepath = 'D:\Data\fepsp_test\lh85_211012_113400';
cd(basepath)

% figure
% stimIdx = 4;
% plot(CattedfepspIO)
% hold on
% yLimit = ylim;
% plot([StimIndexIO{stimIdx}; StimIndexIO{stimIdx}], yLimit, '--k')

% step 1
traces = fepsp_org2traces('data_samples', CattedfepspIO, 'fs', fs,...
    'protocol_id', protocol_id, 'stim_locs', StimIndexIO);

figure, plot(traces{stimIdx})

% step 2
marking_win = fepsp_markings("traces", traces, "fs", fs,...
    "protocol_id", protocol_id,"base_folder", basepath,...
    "intens", intens, "traces_Xlimit", [], "traces_Ylimit", [],...
    "dt", 2, "max_time_tol", 3, "fast_mark", true);





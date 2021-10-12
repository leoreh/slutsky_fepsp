function [AllWindows] = GetWindowsFromProtocol(Protocol,fs,dt,DispWindow)
% [AllWindows] = GetWindowsFromProtocol(Protocol,fs,dt,DispWindow)
% Create Windows by Protocol.
%   Does Not validate inputs! All input are Required.
%   Inputs:
%       Protocol:
%           Required. Stimulation Protocol. Struct From GetProtocol. see
%           GetProtocol for more info.
%       fs:
%           Required. Sampling Frequency, in Hz. Must be a nomeric scalar.
%       dt:
%           Required. Numeric scalar. Dead time between stimuli & earliest
%           possiable responce, in millisec - depends on the synapse you
%           are recourding. Assuming there is nothing intresting dt time
%           before every stimuli and after - therfore AnalyseWin won't
%           include that time. AnalyseWin is used only for initial position
%           of GUI lines and setting ylimits - therefore it does not have
%           much effect.
%       DispWindow:
%           Required, numeric 2 elem. Time in millisec to display
%           in th GUI (xlim). if empty, will take from protocol. In order
%           to be a helpful window, end should be after last responc is
%           expected and start before first stim.
%   Output:
%       AllWindows:
%           Struct with fields:
%               AnalyseWin - row vec, contain IND of all samples in
%                   intresting area (where the responce is expected).
%               AnalyseWinBase - Protocol.nstimX2 numeric, contain starts
%                   and ends of every part of AnalyseWin. They are not
%                   continuous in order to drop stim artifact & dead time. 
%                   If Protocol.nstim > 1 - calcualted as:
%   [1StimStart:PeriodStims:LastStimStart+dt, 1StimEnd:PeriodStims:LastStimEnd-dt]
%                   else
% 	[1StimStart+dt (StimLatency*2)+LastStimEnd-dt]
%           (For Protocol.nstim == 1 case end of window is kind of
%           arbitrary, simpliy going for ~1/3 before & 2/3 after stimulos
%
%               DispWindow - row vec, contain IND of all samples in
%                   display area
%               DispWindowBase - row vec 2 elem numeric, contain starts
%                   and end of display area, same as input if give, else:
%                   If Protocol.nstim > 1 - calcualted as:
%                   [1StimStart-dt LastStimEnd+dt]
%                   else
%                   [1StimStart-dt (StimLatency*2)+LastStimEnd+dt]
%
%               BaselineWindow - row vec, contain IND of all samples in
%                   baseline area (before 1st stim, outside of display) -
%                   between:
%                   -Protocol.StimLatency & Protocol.StimLatency-dt 

% Does Not Get OE STP fs error!
if Protocol.nstim > 1
    PrStims = (1*1000)/Protocol.frStim + Protocol.StimLength/1000; %[mS]
    LastStimStart = (Protocol.nstim-1)*PrStims+Protocol.StimLatency; %[mS]
    
    ANWindowTime = [Protocol.StimLatency:PrStims:LastStimStart;...
        (Protocol.StimLatency+PrStims):PrStims:(LastStimStart+PrStims)]';
    %[1StimStart:PeriodStims:LastStimStart, 1StimEnd:PeriodStims:LastStimEnd]
%     ANWindowTime(end,:) = []; %Created 1 extra stim - remove it
else
    ANWindowTime = [Protocol.StimLatency (Protocol.StimLatency*2)+Protocol.StimLength/1000];
    % arbitrary end limit - Seem comfortable, 1/3 before, two thirds after
end  
AllWindows.AnalyseWinBase = round([ANWindowTime(:,1)+dt,ANWindowTime(:,2)-dt] * fs / 1000);
AllWindows.AnalyseWin = [];
for ii = 1:Protocol.nstim
    AllWindows.AnalyseWin = [AllWindows.AnalyseWin, AllWindows.AnalyseWinBase(ii,1):AllWindows.AnalyseWinBase(ii,2)];
end
if isempty(DispWindow)
    AllWindows.DispWindowBase = round([ANWindowTime(1,1)-dt, ANWindowTime(end,2)+dt] * fs / 1000);
else
    AllWindows.DispWindowBase = DispWindow;
end
AllWindows.DispWindow = AllWindows.DispWindowBase(1):AllWindows.DispWindowBase(2);
AllWindows.BaselineWindow = 1:round((Protocol.StimLatency-dt)* fs / 1000);
end
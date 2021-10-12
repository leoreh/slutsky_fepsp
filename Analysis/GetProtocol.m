function Protocol = GetProtocol(Protocol)
% Protocol = GetProtocol(Protocol)
% Validate & Create Protocol.
%   Inputs:
%       Protocol: Required.
%           Option 1: Text scalar - io/stp, to use our default protocols.
%         To use custom protocols:
%           Option 2: cell array with 6 elements:
%               {nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency}
%               In Order, no name - value.
%           Option 3: struct, same as output.
%       
%        Notice each field/cell description & units: (all are numeric scalars]:
%             nstim = [#] number of stimuli given in each trace
%             frStim = [Hz] Freaquency of stimuli given in each trace
%             StimLength = [uS] How long is each stimuli
%             ProtocolReapeatRate = [Hz] 1 trace = 1 protocol. in what
%                   frequency protocols are reproduced in a long recording,
%                   mesuared from last stimuli given to next trace recording
%                   start.
%             TraceRecordingLenght = [mS] %Total Trace lenght, Including StimLatency & StimLength
%             StimLatency = [mS] time recorded in the trace before the first stimulos is given 
%   Output:
%       Protocol - struct with all fields menthioned above.
%       
%   Examples:
%       Protocol = GetProtocol('io')
%           Use default nstim = 1 protocol.
%       Protocol = GetProtocol({nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency})
%           Validate & Convert to struct
%       Protocol = GetProtocol(ProtocolStruct)
%           Validate Protocol fields
%
% Note - It does not Validate Protocol make sence. Just that fields are
%   numericly okay and are all present.
if ischar(Protocol) || isStringScalar(Protocol)
    Protocol = validatestring(Protocol,{'io','stp','custom'},'fEPSP_analysis_SA','Protocol');
    switch Protocol
        case 'io'
            Protocol = struct();
            Protocol.nstim = 1; %[#]
            Protocol.frStim = inf; %[Hz]
            Protocol.StimLength = 500; %[uS]
%             Protocol.ProtocolReapeatRate = 1/15; %[Hz] %From after Stimuli end
            %Protocol.TraceRecordingLenght = 150; %[mS] %Total, Including StimLatency & StimLength
            Protocol.TraceRecordingLenght = 148.47; % Wierd WCP behavior?
            Protocol.StimLatency = 30; %[mS]
        case 'stp'
            Protocol = struct();
            Protocol.nstim = 3; %[#]
            Protocol.frStim = 50; %[Hz] %From start of stim to start of next one, i.e. disregards stim length
            Protocol.StimLength = 500; %[uS]
%             Protocol.ProtocolReapeatRate = 1/30; %[Hz]
            Protocol.TraceRecordingLenght = 200; %[mS]
            Protocol.StimLatency = 10; %[mS]
        case 'custom'
            AllPrompets = {'Enter number of stimuli given (numeric posative interger)',...
                'Enter the frequency of stimuli given (Hz, numeric posative scalar, 1/peridode between each stim in protocol, if only 1 place ''inf'')',...
                'Enter the length of each stimuli (uS, numeric positive scalar)',...
                'Enter the total length of the recourding (mS, Start to finish including latency)',...
                'Enter the stim latency (mS, How much time passed between start of recourding to first stim given)'};
            UserInputs = str2double(inputdlg(AllPrompets,'Make Custom Protocol'));
            Protocol = struct();
            Protocol.nstim = UserInputs(1); %[#]
            Protocol.frStim = UserInputs(2); %[Hz] %From start of stim to start of next one, i.e. disregards stim length
            Protocol.StimLength = UserInputs(3); %[uS]
            Protocol.TraceRecordingLenght = UserInputs(4); %[mS]
            Protocol.StimLatency = UserInputs(5); %[mS]
    end
elseif iscell(Protocol) && numel(Protocol) == 6
    Tptcl.nstim = Protocol{1}; %[#]
    Tptcl.frStim = Protocol{2}; %[Hz]
    Tptcl.StimLength = Protocol{3}; %[uS]
%     Tptcl.ProtocolReapeatRate = Protocol{4}; %[Hz]
    Tptcl.TraceRecordingLenght = Protocol{5}; %[mS]
    Tptcl.StimLatency = Protocol{6}; %[mS]
    Protocol = Tptcl;
elseif isstruct(Protocol) && isscalar(Protocol)
    if ~all(isfield(Protocol,{'nstim','frStim','StimLength','TraceRecordingLenght','StimLatency'}))
        error(['Protocol must be: char - ''io'' or ''stp'' to use our defaults (see help for info), cell array with 6 elements:'...
            '{nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency}'...
            ', or struct scalar with {nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency} as fields'])
    end
else
    error(['Protocol must be: char - ''io'' or ''stp'' to use our defaults (see help for info), cell array with 6 elements:'...
            '{nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency}'...
            ', or struct scalar with {nstim,frStim,StimLength,ProtocolReapeatRate,TraceRecordingLenght,StimLatency} as fields'])
end

if mod(Protocol.nstim,1)~=0 || Protocol.nstim<=0 || ~isscalar(Protocol.nstim)
    error('When using custom protocol, <nstim> must be a positive interger  scalar [#]')
elseif ~isnumeric(Protocol.frStim) || ~isscalar(Protocol.frStim) || Protocol.frStim<=0 || ~isreal(Protocol.frStim)
    error('When using custom protocol, <frStim> must be a positive scalar [Hz]')
elseif ~isnumeric(Protocol.StimLength) || ~isscalar(Protocol.StimLength) || Protocol.StimLength<=0
    error('When using custom protocol, <StimLength> must be a positive scalar [uS]')
% elseif ~isnumeric(Protocol.ProtocolReapeatRate) || ~isscalar(Protocol.ProtocolReapeatRate) || Protocol.ProtocolReapeatRate<=0
%     error('When using custom protocol, <ProtocolReapeatRate> must be a positive scalar [Hz]')
elseif ~isnumeric(Protocol.TraceRecordingLenght) || ~isscalar(Protocol.TraceRecordingLenght) || Protocol.TraceRecordingLenght<=0
    error('When using custom protocol, <TraceRecordingLenght> must be a positive scalar [mS]')
elseif ~isnumeric(Protocol.StimLatency) || ~isscalar(Protocol.StimLatency) || Protocol.StimLatency<=0
    error('When using custom protocol, <StimLatency> must be a positive scalar [mS]')
end
end

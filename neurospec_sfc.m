function [f, t, c1, sc] = neurospec_sfc(ch, ts, trialType, eventList, eventWin, stepSize, segLength)
%
% usage: 
%
% INPUTS:
%   ch - single element of a channel structure
%   ts - spike timestamps
%   trialType - trial type as defined in getTrialEventParams
%   eventList - list of events within the trials.timestamps structure to be
%       analyzed
%   eventWin - 2-element vector giving the time before and after each event
%       to analyze (i.e., [-1,1] would analyze one second before and after each
%       event)
%   stepSize - size (in seconds) of incremental steps across the event
%       window (for example, 0.1 means the sliding analysis window will advance
%       0.1 s at each step)
%   segLength - size of the analysis window in seconds (for example,
%       stepSize = 0.1 and segLength = 0.2 would mean 200 ms analysis windows
%       get slid across in 100 ms increments)
%   
% VARARGS:
%
% function to calculate spike-field coherence based on the neurospec
% toolbox.
%
% the function sp2a_m1 will be called:
%   function [f,t,cl,sc] = sp2a_m1(sp_type,sp1,dat2,varargin);
%       sp_type = 0, 1, or 2. We will do a type 2 analysis (averaging over
%           repeated trials)
%       sp1 - the spike train, specified as integer numbers of sampling
%           intervals in ascending order
%       dat2 - time series vector
%   varargins (required):
%       1) trig_times - list of trigger times defininig start of each data segment (in samples)
%       2) offset - list of offset values from trigger times to start of
%           each data segment. This should be useful to create sliding
%           windows
%       3) seg_pots - length of each segment, in number of samples
%       4) samp_rate - sampling rate (Hz)
%       5) seg_pwr - segment length in powers of 2
%   varargins (optional):

lfp_root  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal LFPs';
numEvents = length(eventList);

trialEventParams = getTrialEventParams(trialType);
trIdx = extractTrials(ch.trials, trialEventParams);
validTrials = ch.trials(trIdx);

implantID = ch.subject;
if length(implantID) > 5
    implantID = strrep(implantID, '-', '');
end

% load the field potential
[lfp, Fs] = getSingleWireLFP( ch, lfp_root );

for iEvent 1 : numEvents
    eventName = eventList{iEvent};
    event_ts = extractEvent_ts(validTrials, eventName);
    
    % WORKING HERE...
    % IS THERE A FUNDAMENTAL PROBLEM WORKING IN FREQUENCY SPACE, WHERE THE
    % SLIDING WINDOWS SHOULD BE ADJUSTED IN SIZE FOR LOW VS HIGH
    % FREQUENCIES? IS DELTA WITHIN A 200 MS WINDOW EVEN MEANINGFUL?
end

end    % end main function

%**************************************************************************
%**************************************************************************

function [volt_lfp, Fs] = getSingleWireLFP( ch, lfp_root )

implantID = ch.subject;
if length(implantID > 5)
    implantID = strrep(implantID, '-', '');
end

% first, load the field potential
lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
cd(lfp_directory);

sessionDate = ch.date;
if length(sessionDate) > 8
    sessionDate = datestr(sessionDate, 'yyyymmdd');
end

lfp_fileinfo = dir(['*_' sessionDate '.hsdf']);
lfp_fileName = lfp_fileinfo.name;

if length(lfp_fileinfo) ~= 1
    disp([num2str(length(lfp_fileinfo)) ' files found for sessions on ' sessionDate]);
    break;
end

% read in the LFP header and data
header = getHSDHeader( lfp_fileName );
numRecordingSites = header.main.num_channels;
lfp_wireNums = zeros(1, numRecordingSites);
for i_site = 1 : numRecordingSites
    lfp_wireNums(i_site) = header.channel(i_site).original_number;
end

Fs = lfpFs( header );
lfpDuration = getHSDlength( 'filename', lfp_fileName );

lfp = readHSD( lfp_fileName, ...
               numRecordingSites, ...
               header.dataOffset, ...
               Fs, ...
               [0, lfpDuration] );
           
repWire = getRepWire( ch );

lfp_idx = find( lfp_wireNums == repWire );

volt_lfp = int2volt(lfp(lfp_idx, :), 'gain', header.channel(lfp_idx).gain);

end    % function [volt_lfp, Fs] = getSingleWireLFP( ch, lfp_root ) 

%**************************************************************************
%**************************************************************************
function valid_ts = extractEvent_ts( trials, eventName )

valid_ts = [];
for iTr = 1 : length(trials)
    
    if isfield(trials(iTr).timestamps, eventName)
        valid_ts = [valid_ts; trials(iTr).timestamps.(eventName)];
    end
            
end

end    % function valid_ts = extractEvent_ts( trials, eventName ) 
            
            
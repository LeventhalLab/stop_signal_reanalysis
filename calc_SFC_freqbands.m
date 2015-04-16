function [sfc, numSpikes, freqList, numEvents] = calc_SFC_freqbands( ch, sp_ts, trialType, eventList, eventWin, spikeWin, stepSize, varargin )
%
% usage: 
%
% INPUTS:
%   ch - single element of a channel structure
%   sp_ts - spike timestamps
%   trialType - trial type as defined in getTrialEventParams
%   eventList - list of events within the trials.timestamps structure to be
%       analyzed
%   eventWin - 2-element vector giving the time before and after each event
%       to analyze (i.e., [-1,1] would analyze one second before and after each
%       event)
%   spikeWin - duration (in seconds) of the bin size within which to look
%       for spikes to generate the sta
%   stepSize - size (in seconds) of incremental steps across the event
%       window (for example, 0.1 means the sliding analysis window will advance
%       0.1 s at each step)
%   segLength - size of the analysis window in seconds (for example,
%       stepSize = 0.1 and segLength = 0.2 would mean 200 ms analysis windows
%       get slid across in 100 ms increments)
%
% OUTPUTS:

% lfp_root  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal LFPs';
hilbert_directory = sprintf('/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins');

rel_sta_win = 3;    % number of periods over which to calculate 
numEventTypes = length(eventList);

numSteps = floor(range(eventWin) / stepSize);   % number of spike windows contained within a single peri-event window

for iarg = 1 : 2 : nargin - 7
    switch lower(varargin{iarg})
%         case 'lfp_root',
%             lfp_root = varargin{iarg + 1};
        case 'stawin',
            rel_sta_win = varargin{iarg + 1};
        case 'hilbert_directory',
            hilbert_directory = varargin{iarg + 1};
    end
end

if length(sp_ts) > size(sp_ts, 1)   % make sure sp_ts is a column vector
    sp_ts = sp_ts';
end

trialEventParams = getTrialEventParams(trialType);
trIdx = extractTrials(ch.trials, trialEventParams);
validTrials = ch.trials(trIdx);

implantID = ch.subject;
if length(implantID) > 5
    implantID = strrep(implantID, '-', '');
end

subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
hilbert_sessionDir = fullfile(subject_hilbertDir, ch.session);

metadata_filename = [ch.session 'hilbert_metadata.mat'];
metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
load(metadata_filename);
freqList = mean(metadata.freqBands, 2);
numFreqs = length(freqList);
Fs = metadata.Fs;

hilbert_name = ['analytic_' ch.name '.bin'];
hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);

sfc = zeros(numEventTypes, numSteps, numFreqs);
numSpikes = zeros(numEventTypes, numSteps);
validSpikes = cell(numEventTypes, numSteps);

for iEventType = 1 : numEventTypes
    eventName = eventList{iEventType};
    event_ts = extractEvent_ts( validTrials, eventName );
    numEvents = length(event_ts);
%         startSamps = round((event_ts + eventWin(1)) * metadata.Fs);
    eventStartTimes = event_ts + eventWin(1);
    
    for iStep = 1 : numSteps
%         tic
        validSpikes{iEventType, iStep} = [];
        winStartTimes = eventStartTimes + (iStep-1) * stepSize;
        winEndTimes   = winStartTimes + spikeWin;
        for i_singleEvent = 1 : numEvents
            validSpikes{iEventType, iStep} = [validSpikes{iEventType, iStep}; sp_ts(sp_ts > winStartTimes(i_singleEvent) & sp_ts < winEndTimes(i_singleEvent))];
        end
%         toc
        numSpikes(iEventType, iStep) = length(validSpikes{iEventType, iStep});

    end
        
end
            
            
for iFreq = 1 : numFreqs
%     iFreq
    % load the analytic signal
    ansig = readAnalyticSignal(hilbert_name, metadata, [0, metadata.duration], iFreq);
    if length(ansig) > size(ansig, 1)
        ansig = ansig';    % make sure ansig is a column vector
    end
%     y.lfp = real(ansig);
%     y.Fs  = metadata.Fs;
    
    freq = freqList(iFreq);
    sta_width_sec = rel_sta_win / freq;   % sta width in seconds
    sta_width_samps = ceil(sta_width_sec * Fs);    % sta width in samples
    sta_samps_pre   = ceil(sta_width_samps/2);

    for iEventType = 1 : numEventTypes
        
%         eventName = eventList{iEventType};
%         event_ts = extractEvent_ts( validTrials, eventName );
%         numEvents = length(event_ts);
% %         startSamps = round((event_ts + eventWin(1)) * metadata.Fs);
%         eventStartTimes = event_ts + eventWin(1);
        
        for iStep = 1 : numSteps
%             tic
%             validSpikes = [];
%             winStartTimes = eventStartTimes + (iStep-1) * stepSize;
%             winEndTimes   = winStartTimes + spikeWin;
%             for i_singleEvent = 1 : numEvents
%                 validSpikes = [validSpikes; sp_ts(sp_ts > winStartTimes(i_singleEvent) & sp_ts < winEndTimes(i_singleEvent))];
%             end
%             toc
            % validSpikes should now contain all spikes within a
            % specific time window around the current behavioral event
            sp_trig_ansig = zeros(sta_width_samps, length(validSpikes{iEventType, iStep}));
            for iSpike = 1 : length(validSpikes{iEventType, iStep})
                startSamp = round(validSpikes{iEventType, iStep}(iSpike) * Fs ) - sta_samps_pre;
                endSamp = startSamp + sta_width_samps - 1;
                sp_trig_ansig(:, iSpike) = (ansig(startSamp : endSamp));
            end
            
            sp_trig_power = sum(abs(sp_trig_ansig).^2, 1) / sta_width_sec;
            sta_power = sum(abs(mean(sp_trig_ansig, 2)).^2, 1) / sta_width_sec;
            sfc(iEventType, iStep, iFreq) = sta_power / mean(sp_trig_power);
            
        end
        
    end
    
end

end    % function sfc = calc_SFC_freqbands( ch, sp_ts, trialType, eventList, eventWin, spikeWin, stepSize, varargin )

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
    return;
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

%**************************************************************************
%**************************************************************************
            
            
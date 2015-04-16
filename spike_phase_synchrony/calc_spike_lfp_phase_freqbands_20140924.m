function [plv, ppc, numSpikes, spikePhases, freqList, numEvents] = calc_spike_lfp_phase_freqbands_20140924( ch, ...
                                                                                                            sp_ts, ...
                                                                                                            trialType, ...
                                                                                                            eventList, ...
                                                                                                            eventWin, ...
                                                                                                            spikeWin, ...
                                                                                                            stepSize, ...
                                                                                                            unitSaveName, ...
                                                                                                            spikeLFP_metadata, ...
                                                                                                            varargin )
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
%   plv - phase locking value. m x n x p matrix, where m is the number of
%       events, n is the list of frequencies, and p is the number of time
%       steps through individual trials
%   ppc - pairwise phase consistency. m x n x p matrix, where m is the number of
%       events, n is the list of frequencies, and p is the number of time
%       steps through individual trials
%   numSpikes - 
%   spikePhases - 
%   freqList - 
%   numEvents - 

% lfp_root  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal LFPs';
hilbert_directory_1Hz = sprintf('/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins');
hilbert_directory_025Hz = sprintf('/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins');

rel_sta_win = 3;    % number of periods over which to calculate 
numEventTypes = length(eventList);

numSteps = floor(range(eventWin) / stepSize);   % number of spike windows contained within a single peri-event window

for iarg = 1 : 2 : nargin - 9
    switch lower(varargin{iarg})
%         case 'lfp_root',
%             lfp_root = varargin{iarg + 1};
        case 'stawin',
            rel_sta_win = varargin{iarg + 1};
        case 'hilbert_directory_1Hz',
            hilbert_directory_1Hz = varargin{iarg + 1};
        case 'hilbert_directory_025Hz',
            hilbert_directory_025Hz = varargin{iarg + 1};
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

subject_hilbertDir_1Hz = fullfile(hilbert_directory_1Hz, [implantID '_hilbert']);
hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, ch.session);
subject_hilbertDir_025Hz = fullfile(hilbert_directory_025Hz, [implantID '_hilbert']);
hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, ch.session);

metadata_filename = [ch.session 'hilbert_metadata.mat'];
metadata_filename_1Hz = fullfile(hilbert_sessionDir_1Hz, metadata_filename);
metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename);
md_1Hz = load(metadata_filename_1Hz);
md_025Hz = load(metadata_filename_025Hz);

centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
centerFreqs_1Hz = mean(md_1Hz.metadata.freqBands, 2);
maxFreq_025Hz = max(centerFreqs_025Hz);
startFreqIdx_1Hz = find((centerFreqs_1Hz > maxFreq_025Hz), 1, 'first');

freqList = [centerFreqs_025Hz; centerFreqs_1Hz(startFreqIdx_1Hz:end)];   %mean(metadata.freqBands, 2);
num_025HzFreqs = length(centerFreqs_025Hz);

numFreqs = length(freqList);
Fs = md_1Hz.metadata.Fs;

hilbert_name = ['analytic_' ch.name '.bin'];
hilbert_name_1Hz = fullfile(hilbert_sessionDir_1Hz, hilbert_name);
hilbert_name_025Hz = fullfile(hilbert_sessionDir_025Hz, hilbert_name);

if exist(unitSaveName, 'file')
    load(unitSaveName);
else
    plv         = zeros(numEventTypes, numSteps, numFreqs);    % phase locking value
    ppc         = zeros(numEventTypes, numSteps, numFreqs);    % pairwise phase consistency
    numSpikes   = zeros(numEventTypes, numSteps);
    validSpikes = cell(numEventTypes, numSteps);

    spikeLFP_metadata.freqList = freqList;
    spikeLFP_metadata.totalSpikes = length(sp_ts);
end

for iEventType = 1 : numEventTypes
    eventName = eventList{iEventType};
    event_ts = extractEvent_ts( validTrials, eventName );
    numEvents = length(event_ts);
    spikeLFP_metadata.numEvents = numEvents;
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

i_1Hz_freqs = 0;
spikePhases = cell(numEventTypes, numFreqs, numSteps);
tic
startFreq = spikeLFP_metadata.numFreqsComplete + 1;
for iFreq = startFreq : numFreqs

    disp(sprintf('%s, %d of %d frequencies (%f Hz)', ...
                 trialType, ...
                 iFreq, ...
                 numFreqs, ...
                 freqList(iFreq)));
    % load the analytic signal, but only around the relevant events
    if iFreq <= num_025HzFreqs
        hname = hilbert_name_025Hz;
        freqIdx = iFreq;
        metadata = md_025Hz.metadata;
    else
        hname = hilbert_name_1Hz;
        metadata = md_1Hz.metadata;
        freqIdx = startFreqIdx_1Hz + i_1Hz_freqs;
        i_1Hz_freqs = i_1Hz_freqs + 1;
    end
    ansig = readAnalyticSignal(hname, metadata, [0, metadata.duration], freqIdx);
    if length(ansig) > size(ansig, 1)
        ansig = ansig';    % make sure ansig is a column vector
    end

    if spikeLFP_metadata.eventTypesComplete == numEventTypes
        spikeLFP_metadata.eventTypesComplete = 0;
    end
    startEventType = spikeLFP_metadata.eventTypesComplete + 1;
    for iEventType = startEventType : numEventTypes

        for iStep = 1 : numSteps
            
            % validSpikes should now contain all spikes within a
            % specific time window around the current behavioral event
            
            % get sample indices for each spike timestamp
            spike_ts  = validSpikes{iEventType, iStep};
            sample_ts = round(spike_ts * Fs);
            spikePhases{iEventType, iFreq, iStep} = angle(ansig(sample_ts));
            
            plv(iEventType, iStep, iFreq) = abs(sum(exp(1i*spikePhases{iEventType, iFreq, iStep}))) / ...
                                            length(spikePhases{iEventType, iFreq, iStep});
                                        
%             ppc(iEventType, iStep, iFreq) = calc_ppc_20140923(spikePhases{iEventType, iFreq, iStep});
            
        end
        
    end    % for iEventType...
    
    spikeLFP_metadata.eventTypesComplete = numEventTypes;
    spikeLFP_metadata.numFreqsComplete = spikeLFP_metadata.numFreqsComplete + 1;
    save(unitSaveName, 'plv', 'numSpikes', 'spikePhases', 'spikeLFP_metadata');
    
    spikeLFP_metadata.eventTypesComplete = 0;
        
end
toc

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
            
            
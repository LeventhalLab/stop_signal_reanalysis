function [mrv, low_freqs, high_freqs] = phase_amp_Canolty_mrl_20140324(ch, ...
                                                                       eventList, ...
                                                                       eventtWin, ...
                                                                       trialType, ...
                                                                       varargin)
                                                                     
% USAGE:
%
% INPUTS:
%   ch - single element of a channel DB structure
%   eventList - string or cell array of strings containing the names of
%       behavioral events around which to calculate phase-amplitude
%       coupling
%   eventtWin - 2-element vector containing the limits of the window within
%       to look around each event (i.e., [-1 1] would look one second
%       before and after each event)
%   analysisWin - width (in seconds) of each analysis window. For exmaple,
%       if eventtWin = [-1 1] and analysisWin = 0.1, then the first
%       analysis window would be from -1 to 0.9 seconds before the event
%   stepSize - how far to step to generate each point. For example, if
%       eventtWin = [-1 1], analysisWin = 0.1, and stepSize = 0.05 then
%       the modulation index (MI) will be computed in 100 ms windows
%       beginning 1 second before the event, with the window advanced 50 ms
%       for each computation
%   trialType - trial type to analyze (e.g., 'correctgo', etc.)
%
% OUTPUTS:
%   mrv - modulation index. m x n x p x q matrix, where m is the number of
%       event types, n is the number of windowing steps around each event,
%       p is the number of "phase" frequencies, and q is the number of
%       "amplitude" frequencies
%   surrogate - 

% algorithm to compute phaes-amplitude coupling based on the algorithm of
% Canolty et al, Science, 2006. This differs from the phase_amp_comodugram
% function in that this version computes the modulation index (MI) for
% discrete windows across trials and averages them together. This version
% can be used to examine only relationships during trials, or in sliding
% windows around discrete behavioral events.

%Can use a shuffle test to find areas of significant MI elevation?

% rev history
% 3/24/2014 - allowed the "low frequencies" to include data read from the
%             1Hz spaced analytic signal files
low_freq_range  = [0, 21];
high_freq_range = [8, 100];

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
lfp_root          = '/Volumes/Recordings/dan/Leventhal Neuron 2012_summary/Leventhal Neuron 2012_LFPs';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';

if ~iscell(eventList)
    eventList = {eventList};
end
for iarg = 1 : 2 : nargin - 4
    switch lower(varargin{iarg})
        case 'lowfreqrange',
            low_freq_range = varargin{iarg + 1};
        case 'highfreqrange',
            high_freq_range = varargin{iarg + 1};
        case 'phasebins',
            phasebins = varargin{iarg + 1};
        case 'hilbert1hzdir',
            hilbert_1Hz_directory = varargin{iarg + 1};
        case 'hilbert025hzdir',
            hilbert_025Hz_directory = varargin{iarg + 1};
    end
end

% navigate to the directory in which the analytic signals are stored
implantID = implantID_from_ratID(ch.name(1:5));
subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);

hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, ch.session);
hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, ch.session);

metadata_filename_1Hz   = [ch.session 'hilbert_metadata.mat'];
metadata_filename_025Hz = [ch.session 'hilbert_metadata.mat'];
metadata_filename_1Hz   = fullfile(hilbert_sessionDir_1Hz, metadata_filename_1Hz);
metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename_025Hz);

if ~exist(metadata_filename_1Hz, 'file')
    error('metadata for 1 Hz bands not found.');
end
if ~exist(metadata_filename_025Hz, 'file')
    error('metadata for 0.25 Hz bands not found.');
end

md_1Hz   = load(metadata_filename_1Hz);
md_025Hz = load(metadata_filename_025Hz);
Fs = md_1Hz.metadata.Fs;

% hilbert_name_1Hz   = ['analytic_' ch.name '.bin'];
% hilbert_name_025Hz = ['analytic_' ch.name '.bin'];
% hilbert_name_1Hz   = fullfile(hilbert_sessionDir_1Hz, hilbert_name_1Hz);
% hilbert_name_025Hz = fullfile(hilbert_sessionDir_025Hz, hilbert_name_025Hz);

centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);

max025Hz = max(centerFreqs_025Hz);

low_freq_idx_025  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
low_freq_idx_1    = find(centerFreqs_1Hz >= max025Hz & centerFreqs_1Hz <= low_freq_range(2));
low_freq_idx      = [low_freq_idx_025; low_freq_idx_1];
high_freq_idx     = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

low_freqs  = [centerFreqs_025Hz(low_freq_idx_025); centerFreqs_1Hz(low_freq_idx_1)];
high_freqs = centerFreqs_1Hz(high_freq_idx);

num_low_freq_025 = length(low_freq_idx_025);
num_low_freq  = length(low_freq_idx);
num_high_freq = length(high_freq_idx);

numEventTypes = length(eventList);
numSamps = round(range(eventtWin) * Fs);

trialEventParams = getTrialEventParams(trialType);
trList = extractTrials(ch.trials, trialEventParams);

% figure out how many VALID trials there are (that is, trials that don't
% occur to close to the beginning or end of the session so that we can't
% extract a full window around all events
numTrials = length(trList);
minValidTrials = numTrials;
for iEventType = 1 : numEventTypes
    eventName = eventList{iEventType};
    numValidTrials = 1;
    for iTr = 1 : numTrials
        if ~isfield(ch.trials(trList(numValidTrials)).timestamps, eventName)
            if numValidTrials == 1
                trList = trList(2:end,:);
            elseif numValidTrials == length(trList)
                trList = trList(1:numValidTrials);
            else
                trList = [trList(1:numValidTrials-1); trList(numValidTrials+1:end)];
            end
            continue
        end
        event_ts = ch.trials(trList(numValidTrials)).timestamps.(eventName);
        
        tlim = event_ts + eventtWin;
        % make sure the time window doesn't extend before the start of the
        % recording or after the end of the recording
        if any(tlim < 0) || any(tlim > md_1Hz.metadata.duration)
            if numValidTrials == 1
                trList = trList(2:end,:);
            elseif numValidTrials == length(trList)
                trList = trList(1:numValidTrials);
            else
                trList = [trList(1:numValidTrials-1); trList(numValidTrials+1:end)];
            end
            continue;
        end
        numValidTrials = numValidTrials + 1;
    end
    numTrials = numValidTrials - 1;
    if numValidTrials - 1 < minValidTrials
        minValidTrials = numValidTrials - 1;
    end
end
numValidTrials = minValidTrials;
        
mrv = zeros(numEventTypes, num_low_freq, num_high_freq, numSamps);
for iEventType = 1 : numEventTypes
    iEventType
    
    as2 = zeros(num_high_freq, numValidTrials, numSamps);
    as1 = zeros(num_low_freq, numValidTrials, numSamps);
    
    for i_f2 = 1 : num_high_freq
% i_f2
%         disp(sprintf('i_f2 = %d', i_f2))
        as2(i_f2, :, :) = getAnalyticAroundEvent_20140916( ch, ...
                                                           trList, ...
                                                           high_freq_idx(i_f2), ...
                                                           eventList{iEventType}, ...
                                                           eventtWin, ...
                                                           'hilbertdir', hilbert_1Hz_directory );
    end
    
    sig_amp = abs(as2);
    
    
    for i_f1 = 1 : num_low_freq
%         disp(sprintf('i_f1 = %d', i_f1))
        if i_f1 <= num_low_freq_025 
            activeHilbertDir = hilbert_025Hz_directory;
        else
            activeHilbertDir = hilbert_1Hz_directory;
        end
        
        as1(i_f1, :, :) = getAnalyticAroundEvent_20140916( ch, ...
                                                           trList, ...
                                                           low_freq_idx(i_f1), ...
                                                           eventList{iEventType}, ...
                                                           eventtWin, ...
                                                           'hilbertdir', activeHilbertDir );
                                  
    end
    
    phase_angles = angle(as1);
    phase_angles(phase_angles < 0) = phase_angles(phase_angles < 0) + 2 * pi;
    
%     surrogate_mean = zeros(num_low_freq, num_high_freq);
%     surrogate_std  = zeros(num_low_freq, num_high_freq);
    for i_f1 = 1 : num_low_freq
        for i_f2 = 1 : num_high_freq
            if low_freqs(i_f1) >= high_freqs(i_f2); continue; end
            
            phaseAmpfunction = squeeze(sig_amp(i_f2,:,:)) .* squeeze(exp(1i*phase_angles(i_f1,:,:)));
            mrv(iEventType, i_f1, i_f2, :) = mean(phaseAmpfunction, 1);
            
        end
    end
    
end

% surrogate.mean = surrogate_mean;
% surrogate.std  = surrogate_std;

end
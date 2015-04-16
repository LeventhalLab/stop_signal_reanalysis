function [MI, low_freqs, high_freqs] = mean_phase_amp_comodugram(ch, ...
                                                                 eventList, ...
                                                                 eventtWin, ...
                                                                 analysisWin, ...
                                                                 stepSize, ...
                                                                 trialType, ...
                                                                 varargin)
                                                              
% USAGE:
%
% INPUTS:
%   ch - single element of a channel DB structure
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
%   MI - modulation index. m x n x p x q matrix, where m is the number of
%       event types, n is the number of windowing steps around each event,
%       p is the number of "phase" frequencies, and q is the number of
%       "amplitude" frequencies

% algorithm to compute phaes-amplitude coupling based on the algorithm of
% Tort et al, PNAS, 2008. This differs from the phase_amp_comodugram
% function in that this version computes the modulation index (MI) for
% discrete windows across trials and averages them together. This version
% can be used to examine only relationships during trials, or in sliding
% windows around discrete behavioral events.

%Can use a shuffle test to find areas of significant MI elevation?

low_freq_range  = [0, 10];
high_freq_range = [8, 100];

phasebins = (10:20:350) * pi / 180;

bitOrder = 'b';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
lfp_root          = '/Volumes/Recordings/dan/Leventhal Neuron 2012_summary/Leventhal Neuron 2012_LFPs';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';

for iarg = 1 : 2 : nargin - 6
    switch lower(varargin{iarg})
        case 'lowfreqrange',
            low_freq_range = varargin{iarg + 1};
        case 'highfreqrange',
            high_freq_range = varargin{iarg + 1};
        case 'phasebins',
            phasebins = varargin{iarg + 1};
        case 'hilbert1Hzdir',
            hilbert_1Hz_directory = varargin{iarg + 1};
        case 'hilbert025Hzdir',
            hilbert_025Hz_directory = varargin{iarg + 1};
    end
end
nbins = length(phasebins);
binLims = 0 : diff(phasebins(1:2)) : 2*pi;
            
% navigate to the directory in which the analytic signals are stored
implantID = implantID_from_ratID(ch.name(1:3));
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

hilbert_name_1Hz   = ['analytic_' ch.name '.bin'];
hilbert_name_025Hz = ['analytic_' ch.name '.bin'];
hilbert_name_1Hz   = fullfile(hilbert_sessionDir_1Hz, hilbert_name_1Hz);
hilbert_name_025Hz = fullfile(hilbert_sessionDir_025Hz, hilbert_name_025Hz);

centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);

low_freq_idx  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

low_freqs  = centerFreqs_025Hz(low_freq_idx);
high_freqs = centerFreqs_1Hz(high_freq_idx);

num_low_freq  = length(low_freq_idx);
num_high_freq = length(high_freq_idx);

trialEventParams = getTrialEventParams( trialType );
trList = extractTrials( ch.trials, trialEventParams );

numEventTypes = length(eventList);
winSamps      = floor(analysisWin * Fs);
winStartTimes = eventtWin(1) : stepSize : (eventtWin(2) - analysisWin);
winStartSamps = ceil((winStartTimes-eventtWin(1)) * Fs);
if winStartSamps(1) == 0; winStartSamps(1) = 1; end

numSteps = length(winStartTimes);    %floor((range(eventtWin) - analysisWin) / stepSize) + 1;
numSamps = round(range(eventtWin) * Fs);

trialEventParams = getTrialEventParams(trialType);
trList = extractTrials(ch.trials, trialEventParams);
numTrials = length(trList);

MI = zeros(numEventTypes, numSteps, num_low_freq, num_high_freq);
for iEventType = 1 : numEventTypes
    iEventType
    
    as2 = zeros(num_high_freq, numTrials, numSamps);
    for i_f2 = 1 : num_high_freq

        as2(i_f2, :, :) = getAnalyticAroundEvent( ch, ...
                                                  high_freq_idx(i_f2), ...
                                                  eventList{iEventType}, ...
                                                  trialType, ...
                                                  eventtWin, ...
                                                  'hilbertdir', hilbert_1Hz_directory );
    end
    sig_amp = abs(as2);

    for i_f1 = 1 : num_low_freq
        
        as1 = getAnalyticAroundEvent( ch, ...
                                      low_freq_idx(i_f1), ...
                                      eventList{iEventType}, ...
                                      trialType, ...
                                      eventtWin, ...
                                      'hilbertdir', hilbert_025Hz_directory );
        phase_angles = angle(as1);
        phase_angles(phase_angles < 0) = phase_angles(phase_angles < 0) + 2 * pi;

        for i_f2 = 1 : num_high_freq
            
            if centerFreqs_025Hz(low_freq_idx(i_f1)) >= centerFreqs_1Hz(high_freq_idx(i_f2)); continue; end
            
            for iStep = 1 : numSteps
                endSamp = winStartSamps(iStep) + winSamps;
                if endSamp > size(as1, 2); endSamp = size(as1, 2); end
                MI(iEventType, iStep, i_f1, i_f2) = meanAmp_per_phase( phase_angles(:, winStartSamps(iStep) : endSamp), ...
                                                                       squeeze(sig_amp(i_f2, :, winStartSamps(iStep) : endSamp)), ...
                                                                       phasebins);
            end

        end    % for i_f2...
        
    end    % for i_f1...
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MI = meanAmp_per_phase( ph, amp, phasebins )

y = zeros(length(phasebins), 1);
binLims = 0 : diff(phasebins(1:2)) : 2*pi;

for iBin = 1 : length(binLims)-1
    y(iBin) = mean(amp(ph >= binLims(iBin) & ph < binLims(iBin + 1)));
end
U = ones(1, length(phasebins)) / length(phasebins);

% CORRECTION FOR POSSIBLE EMPTY BINS
if any(isnan(y))
    y(isnan(y)) = nanmean(y);
end
y = y ./ nansum(y);    % normalize over sum of mean amplitudes across bins
MI = KLdistance(y, U) / log10(length(phasebins));
        
end
function ansig = getAnalyticAroundEvent_20140916( ch, trList, freqIdx, eventName, twin, varargin )
%
% function to extract the analytic signals around a specific behavioral
% event for a specific trial type for a single "channel" (a
% tetrode-session)
%
% usage: ansig = getAnalyticAroundEvent( ch, freqIdx, eventName, trialType, twin )
%
% INPUTS:
%   ch - a single element of a channel structure
%   freqIdx - index of the frequency of interest in the metadata file
%   eventName - name (case-sensitive) of the behavioral event around which
%       to extract the analytic signals
%   trialType - type of trials to use (see "getTrialEventParams" for a list
%       of possible trial types
%   twin - time window to extract, a 2-element vector (ie, [-1 1] will
%       extract one second around the behavioral event
%
% VARARGINs:
%   'hilbertdir' - directory in which to find the hilbert transformed data;
%       default is '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins'
%
% OUTPUTS:
%   ansig - m x n matrix containing the analytic signals (n samples for
%       each behavioral event)

hilbert_directory = '\\170.20.138.142/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

for iarg = 1 : 2 : nargin - 5
    switch lower(varargin{iarg})
        case 'hilbertdir',
            hilbert_directory = varargin{iarg + 1};
    end
end

implantID = implantID_from_ratID(ch.name(1:5));
subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);

hilbert_sessionDir = fullfile(subject_hilbertDir, ch.session);
metadata_filename = [ch.session 'hilbert_metadata.mat'];
metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
if ~exist(metadata_filename, 'file')
    error([metadata_filename ' could not be found.']);
end
load(metadata_filename);

hilbert_name = ['analytic_' ch.name '.bin'];
hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
if ~exist(hilbert_name, 'file')
    error([hilbert_name ' could not be found.']);
end
fullsig = readAnalyticSignal(hilbert_name, metadata, [0, metadata.duration], freqIdx);
% trialEventParams = getTrialEventParams(trialType);
% trIdx = extractTrials(ch.trials, trialEventParams);

numTrials = length(trList);

samps_per_window = round(range(twin) * metadata.Fs);
ansig = zeros(numTrials, samps_per_window);
ansig = complex(ansig, ansig);

% numValidTrials = 0;
for iTr = 1 : numTrials
%     iTr
    if ~isfield(ch.trials(trList(iTr)).timestamps, eventName); continue; end
    event_ts = ch.trials(trList(iTr)).timestamps.(eventName);
    
    tlim = event_ts + twin;
    % make sure the time window doesn't extend before the start of the
    % recording or after the end of the recording
%     if any(tlim < 0) || any(tlim > metadata.duration)
%         if numValidTrials == 0
%             ansig = ansig(2:end,:);
%         elseif numValidTrials == size(ansig, 1) - 1
%             ansig = ansig(1:numValidTrials, :);
%         else
%             ansig = [ansig(1:numValidTrials,:); ansig(numValidTrials + 1:end, :)];
%         end
%         continue;
%     end
%     numValidTrials = numValidTrials + 1;
    startSamp = round(tlim(1) * metadata.Fs);
    ansig(iTr, :) = fullsig( startSamp : startSamp+samps_per_window-1 );
end
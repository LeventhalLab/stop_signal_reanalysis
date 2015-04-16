function ansig = getAnalyticAround_ts( ch, freqIdx, ts, twin )
%
% function to extract the analytic signals around a set of timesatmps
% for a single "channel" (a tetrode-session). This differs from
% getAnalyticAroundEvent because the timestamps are specified as input to
% this function, rather than pulled out for each event
%
% usage: ansig = getAnalyticAround_ts( ch, freqIdx, ts, twin )
%
% INPUTS:
%   ch - a single element of a channel structure
%   freqIdx - index of the frequency of interest in the metadata file
%   ts - vector of timestamps around which to extract the analytic signal
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
%       each of m behavioral events)

hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

for iarg = 1 : 2 : nargin - 5
    switch lower(varargin{iarg})
        case 'hilbertdir',
            hilbert_directory = varargin{iarg + 1};
    end
end

implantID = implantID_from_ratID(ch.name(1:3));
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

num_ts = length(ts);

samps_per_window = round(range(twin) * metadata.Fs);
ansig = zeros(num_ts, samps_per_window);
ansig = complex(ansig, ansig);

numValid_ts = 0;
for i_ts = 1 : num_ts
    
    tlim = ts(i_ts) + twin;
    % make sure the time window doesn't extend before the start of the
    % recording or after the end of the recording
    if any(tlim < 0) || any(tlim > metadata.duration)
        if numValid_ts == 0
            ansig = ansig(2:end,:);
        elseif numValid_ts == size(ansig, 1) - 1
            ansig = ansig(1:numValid_ts, :);
        else
            ansig = [ansig(1:numValid_ts,:); ansig(numValid_ts + 1:end, :)];
        end
        continue;
    end
    numValid_ts = numValid_ts + 1;
    
    ansig(numValid_ts, :) = readAnalyticSignal(hilbert_name, metadata, tlim, freqIdx);
    
end
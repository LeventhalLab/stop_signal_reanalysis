function createSurrogate_phase_RT_dist( ch, varargin )
%
% function to create the surrogate distribution of phase as a function of
% RT
%
% INPUTS:
%   ch - a single channel of a channel DB structure
%   freq - center frequency of band from which to extract phases
%   hilbert_directory - parent directory for where the hilbert-transformed
%       data are stored
%
% VARARGINs:
%   'iterations' - number of surrogate iterations to calculate
%   'hilbertdirectory' - root directory for analytic signals
%   'phasertcorrdir - root directory for phase-RT correlation calculations

% adapted from VanRullen et al, Ongoing
% EEG phase..., Frontiers in Psychology, 2011; details in Drewes and
% Vanrullen, "This is the rhtyhm of your eyes: the phase of ongoing
% electroencephalogram oscillations modulates saccadic reaction time.", J
% Neurosci, 2011

numIterations = 100;
trialFract    = 0.2;   % fraction of total number of trials to include in each surrogate distribution

hilbert_directory     = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

for iarg = 1 : 2 : nargin - 2
    switch lower(varargin{iarg})
        case 'iterations',
            numIterations = varargin{iarg + 1};
        case 'hilbertdirectory',
            hilbert_directory = varargin{iarg + 1};
        case 'phasertcorrdir',
            phaseRTcorr_directory = varargin{iarg + 1};
    end
end

implantID = implantID_from_ratID(ch.name(1:3));

subject_phaseRTcorrdir = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr']);
phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorrdir, ch.session);
surrogateName = ['phase_RT_surrogates_' ch.name '.mat'];
surrogateName = fullfile(phaseRTcorr_sessionDir, surrogateName);
surrogateMetadataName = ['phase_RT_surrogate_metadata_' ch.name '.mat'];
surrogateMetadataName = fullfile(phaseRTcorr_sessionDir, surrogateMetadataName);
surrogate_metadata_loaded = false;
if exist(surrogateMetadataName, 'file')
    load(surrogateMetadataName);
    surrogate_metadata_loaded = true;
    if surrogate_metadata.completeFreqs == length(surrogate_metadata.freqList)
        disp(['surrogate phase-RT relationships already calculated for ' ch.name]);
        return;
    end
end
if exist(surrogateName, 'file')
    load(surrogateName);
end

phase_RTcorr_name = ['phase_RT_analysis_' ch.name '.mat'];
phase_RTcorr_name = fullfile(phaseRTcorr_sessionDir, phase_RTcorr_name);
if ~exist(phase_RTcorr_name, 'file')
    disp([phase_RTcorr_name ' not found.'])
    return
end

subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
hilbert_sessionDir = fullfile(subject_hilbertDir, ch.session);
hilbert_name       = ['analytic_' ch.name '.bin'];
hilbert_name       = fullfile(hilbert_sessionDir, hilbert_name);
if ~exist(hilbert_name, 'file')
    disp([hilbert_name ' not found.'])
    return
end
hilbert_metadata_name = [ch.session 'hilbert_metadata.mat'];
hilbert_metadata_name = fullfile(hilbert_sessionDir, hilbert_metadata_name);
load(hilbert_metadata_name);
Fs = metadata.Fs;
freqList = mean(metadata.freqBands, 2);
numFreqs = length(freqList);

load(phase_RTcorr_name);
% phaseRTcorr_metadata fields as established in script_RTcorrelations:
%         phaseRTcorr_metadata.freqs = centerFreqs;
%         phaseRTcorr_metadata.chNames = metadata.chNames;
%         phaseRTcorr_metadata.eventList = eventList;
%         phaseRTcorr_metadata.twin = twin;
%         phaseRTcorr_metadata.Fs = metadata.Fs;
%         phaseRTcorr_metadata.trialType = trialType;
%         phaseRTcorr_metadata.RTquantiles = RTquantiles;
%         phaseRTcorr_metadata.RTcutoffs = RTquantile_borders;

eventList = phaseRTcorr_metadata.eventList;
trialType = phaseRTcorr_metadata.trialType;
twin      = phaseRTcorr_metadata.twin;
numEvents = length(eventList);

trialEventParams = getTrialEventParams(trialType);
validTrials      = extractTrials2(ch.trials, trialEventParams);
numTrials        = length(validTrials);

trials_per_iteration = ceil(trialFract * numTrials);

samps_per_window = round(range(twin) * Fs);
mrl = zeros(numFreqs, numEvents, numIterations, samps_per_window);

if ~surrogate_metadata_loaded
    surrogate_metadata.channel = ch.name;
    surrogate_metadata.Fs = Fs;
    surrogate_metadata.trialType = trialType;
    surrogate_metadata.centerFreqs = freqList;
    surrogate_metadata.twin = twin;
    surrogate_metadata.numIterations = numIterations;
    surrogate_metadata.eventList = eventList;
    surrogate_metadata.completeFreqs = 0;
    surrogate_metadata.freqList = freqList;
end

for iFreq = surrogate_metadata.completeFreqs + 1 : numFreqs
    tic
    for iIteration = 1 : numIterations
        
        trialSampleIdx = randperm(numTrials, trials_per_iteration);

        % cycle through events, extract phase angles at each time point to
        % create a surrogate distribution of ITCs
        for iEvent = 1 : numEvents
            
            ts = zeros(1, trials_per_iteration);
            for iTrial = 1 : trials_per_iteration
                ts(iTrial) = validTrials(trialSampleIdx(iTrial)).timestamps.(eventList{iEvent});
            end
            ansig = getAnalyticAround_ts( ch, ...
                                          iFreq, ...
                                          ts, ...
                                          twin );
            freqPhase = angle(ansig);

            mrl(iFreq, iEvent, iIteration, :) = calcMRV(freqPhase);
            
        end

    end
    surrogate_metadata.completeFreqs = surrogate_metadata.completeFreqs + 1;
    save(surrogateName, 'mrl');
    save(surrogateMetadataName, 'surrogate_metadata');
    toc
end


        
function createSurrogate_phase_RT_dist_fail_succ( ch, varargin )
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
surrogateName = ['phase_RT_surrogates_fail_succ' ch.name '.mat'];
surrogateName = fullfile(phaseRTcorr_sessionDir, surrogateName);
surrogateMetadataName = ['phase_RT_surrogate_metadata_fail_succ' ch.name '.mat'];
surrogateMetadataName = fullfile(phaseRTcorr_sessionDir, surrogateMetadataName);
surrogate_metadata_loaded = false;
if exist(surrogateMetadataName, 'file')
    load(surrogateMetadataName);
    surrogate_metadata_loaded = true;
    if surrogate_metadata_fail_succ.completeFreqs == length(surrogate_metadata_fail_succ.freqList)
        disp(['surrogate phase-RT relationships already calculated for ' ch.name]);
        return;
    end
end
if exist(surrogateName, 'file')
    load(surrogateName);
end

phase_RTcorr_name = ['phase_RTanalysisStopSuccessFail_' ch.name '.mat'];
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
%         phaseRTcorrStopSuccessFail_metadata.freqs = centerFreqs;
%         phaseRTcorrStopSuccessFail.chNames = metadata.chNames;
%         phaseRTcorrStopSuccessFail.eventList = eventList;
%         phaseRTcorrStopSuccessFail.twin = twin;
%         phaseRTcorrStopSuccessFail.Fs = metadata.Fs;
%         phaseRTcorrStopSuccessFail.trialTypeSucc = trialType;
%         phaseRTcorrStopSuccessFail.trialTypeFail = trialType2;
%         phaseRTcorrStopSuccessFail.RTquantiles = RTquantiles;
%         phaseRTcorrStopSuccessFail.RTcutoffs = RTquantile_borders;

eventList = phaseRTcorrStopSuccessFail_metadata.eventList;
trialTypeSucc = phaseRTcorrStopSuccessFail_metadata.trialTypeSucc;
trialTypeFail = phaseRTcorrStopSuccessFail_metadata.trialTypeFail;
twin      = phaseRTcorrStopSuccessFail_metadata.twin;
numEvents = length(eventList);

trialEventParamsSucc = getTrialEventParams(trialTypeSucc);
trialEventParamsFail = getTrialEventParams(trialTypeFail);
validTrialsSucc      = extractTrials2(ch.trials, trialEventParamsSucc);
validTrialsFail      = extractTrials2(ch.trials, trialEventParamsFail);
validTrials = [validTrialsSucc;validTrialsFail];
numTrials        = length(validTrials);

trials_per_iteration = ceil(trialFract * numTrials);

samps_per_window = round(range(twin) * Fs);
mrl = zeros(numFreqs, numEvents, numIterations, samps_per_window);

if ~surrogate_metadata_loaded
    surrogate_metadata_fail_succ.channel = ch.name;
    surrogate_metadata_fail_succ.Fs = Fs;
    surrogate_metadata_fail_succ.trialTypeSucc = trialTypeSucc;
    surrogate_metadata_fail_succ.trialTypeFail = trialTypeFail;
    surrogate_metadata_fail_succ.centerFreqs = freqList;
    surrogate_metadata_fail_succ.twin = twin;
    surrogate_metadata_fail_succ.numIterations = numIterations;
    surrogate_metadata_fail_succ.eventList = eventList;
    surrogate_metadata_fail_succ.completeFreqs = 0;
    surrogate_metadata_fail_succ.freqList = freqList;
end

for iFreq = surrogate_metadata_fail_succ.completeFreqs + 1 : numFreqs
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
    surrogate_metadata_fail_succ.completeFreqs = surrogate_metadata_fail_succ.completeFreqs + 1;
    save(surrogateName, 'mrl');
    save(surrogateMetadataName, 'surrogate_metadata_fail_succ');
    toc
end


        
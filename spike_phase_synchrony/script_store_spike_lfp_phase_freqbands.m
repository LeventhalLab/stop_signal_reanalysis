% script_store_spike_lfp_phase_freqbands

% script calls calc_spike_lfp_phase_freqbands_20140207, which calculates
% the phase locking value (plv) and pairwise phase coherence (ppc). For
% references, see:

% van Wingerden et al, "Learning-Associated Gamma-Band Phase-Locking of
%   Action-Outcome Selective Neurons in Orbitofrontal Cortex", J Neurosci,
%   2010
% Vinck et al, "The pairwise phase consistency: A bias-free
%   measure of rhythmic neuronal synchronization", Neuroimage, 2010

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
spikeLFP_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/spikeLFP_windowed';
[chDB_list, chDB_fnames, spikeDB_list, spikeDB_fnames] = get_chStructs_for_analysis;

eventList{1} = {'noseCenterIn'};
eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
eventList{3} = eventList{1};
eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

numTrialTypes = length(trialTypeList);
eventtWin   = zeros(numTrialTypes, 2);
analysisWin = zeros(numTrialTypes, 1);
stepSize    = zeros(numTrialTypes, 1);
eventtWin(1,:) = [-1 2];   % for analysis of all trials
analysisWin(1) = 3;
stepSize(1)    = 3;
for iTrialType = 2 : numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
    analysisWin(iTrialType) = 0.1;
    stepSize(iTrialType)    = 0.05;
end

for i_chDB = 1 : 1%length(chDB_list)
    
    % first, load the relevant channel and spike DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    if ~exist(spikeDB_list{i_chDB}, 'var')
        spikeDB_file = fullfile(chDB_directory, spikeDB_fnames{i_chDB});
        disp(['loading ' spikeDB_file]);
        load( spikeDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir_1Hz, 'dir')
        disp([subject_hilbertDir_1Hz ' not found. Skipping ' implantID '...'])
        continue
    end
    if ~exist(subject_hilbertDir_025Hz, 'dir')
        disp([subject_hilbertDir_025Hz ' not found. Skipping ' implantID '...'])
        continue
    end

    subject_spikeLFPdir = fullfile(spikeLFP_directory, [implantID '_spike_phase']);
    if ~exist(subject_spikeLFPdir, 'dir')
        mkdir(subject_spikeLFPdir);
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
        spikeDB_info = whos( [spikeDB_list{i_chDB}(1:3) 'spike'] );
    else
        implantID = chDB_list{i_chDB}(1:5);
        chDB_info = whos( [implantID 'Ch*'] );
        spikeDB_info = whos( implantID );
    end
    channels = eval( chDB_info.name );
    spike    = eval( spikeDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iTrialType = 1 : length(trialTypeList)
        trialType = trialTypeList{iTrialType};

        spikeLFP_metadata.eventList           = eventList{iTrialType};
        spikeLFP_metadata.trialType           = trialType;
        spikeLFP_metadata.eventtWin           = eventtWin(iTrialType, :);
        spikeLFP_metadata.analysisWin         = analysisWin(iTrialType);
        spikeLFP_metadata.stepSize            = stepSize(iTrialType);
        spikeLFP_metadata.eventTypesComplete  = 0;
        spikeLFP_metadata.numFreqsComplete    = 0;
        
        numEventTypes = length(eventList{iTrialType});
        
        for iSession = 1 : numSessions
            
            disp(sprintf('%s, %s', trialType, sessionList{iSession}));
            
            hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
            hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
            if ~exist(hilbert_sessionDir_1Hz, 'dir')
                disp([hilbert_sessionDir_1Hz ' not found. Skipping ' sessionList{iSession} '...'])
                continue
            end
            if ~exist(hilbert_sessionDir_025Hz, 'dir')
                disp([hilbert_sessionDir_025Hz ' not found. Skipping ' sessionList{iSession} '...'])
                continue
            end
            
            session_spikeLFPdir = fullfile(subject_spikeLFPdir, sessionList{iSession});
            if ~exist(session_spikeLFPdir, 'dir')
                mkdir(session_spikeLFPdir);
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            session_spikeList = extractChannels(cp, spike);
            sessionSpikes     = spike(session_spikeList);
            
%             session_chList = extractChannels( cp, channels );
%             sessionChannels = channels(session_chList);
            numSessionUnits = length(sessionSpikes);
            
            for iUnit = 1 : numSessionUnits
                
                disp(sprintf('unit %d of %d', iUnit, numSessionUnits));
                
                sp = sessionSpikes{iUnit};
%                 spName{iUnit} = sp.name;
                
                unitSaveName = ['spikeLFP_phase_' sp.name '_' trialType '.mat'];
                unitSaveName = fullfile(session_spikeLFPdir, unitSaveName);

                cp = initChanParams();
                cp.tetrode = sp.tetrode.name;
                cp.session = sp.session;

                chList = extractChannels(cp, channels);
                if isempty(chList); continue; end

                ch = channels{chList};
                sp_ts = sp.timestamps.spike;
                
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

                if exist(unitSaveName, 'file')
                    load(unitSaveName);
                    if (spikeLFP_metadata.eventTypesComplete == numEventTypes) && ...
                       (spikeLFP_metadata.numFreqsComplete == numFreqs)
                        continue
                    end
                end
                
                [plv, ppc, numSpikes, spikePhases, freqList, numEvents] = calc_spike_lfp_phase_freqbands_20140924( ch, ...
                                                                                                                   sp_ts, ...
                                                                                                                   trialType, ...
                                                                                                                   eventList{iTrialType}, ...
                                                                                                                   eventtWin(iTrialType, :), ...
                                                                                                                   analysisWin(iTrialType), ...
                                                                                                                   stepSize(iTrialType), ...
                                                                                                                   unitSaveName, ...
                                                                                                                   spikeLFP_metadata, ...
                                                                                                                   'hilbert_directory_1Hz', hilbert_1Hz_directory, ...
                                                                                                                   'hilbert_directory_025Hz', hilbert_025Hz_directory );
                
                
            end    % for iUnit...
            
        end    % for iSession...
        
    end    % for iTrialType...
    
end    % for i_chDB
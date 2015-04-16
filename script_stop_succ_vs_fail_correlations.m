% script_RTcorrelations

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};
numEvents = length(eventList);
twin = [-1 1];

% parameters for phase-RT analysis; adapted from VanRullen et al, Ongoing
% EEG phase..., Frontiers in Psychology, 2011; details in Drewes and
% Vanrullen, "This is the rhtyhm of your eyes: the phase of ongoing
% electroencephalogram oscillations modulates saccadic reaction time.", J
% Neurosci, 2011
numQuantiles = 2;
trialType = 'correctstop';
trialType2 = 'failedstop';


for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir, 'dir')
        disp([subject_hilbertDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    subject_powerRTcorrdir = fullfile(powerRTcorr_directory, [implantID '_powerRTcorrStopSuccessFail']);
    if ~exist(subject_powerRTcorrdir, 'dir')
        mkdir(subject_powerRTcorrdir);
    end
    subject_phaseRTcorrdir = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorrStopSuccessFail']);
    if ~exist(subject_phaseRTcorrdir, 'dir')
        mkdir(subject_phaseRTcorrdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    [RT, ~, sessionList] = collect_RT_MT_by_rat(channels, trialType);
    [RT2, ~, ~] = collect_RT_MT_by_rat(channels, trialType2);
    numSessions = length( sessionList );
    
    
    
    powerRTcorrStopSuccessFail = cell(1, numSessions);    % cell array to store correlations between continuous narrow-band power and RT
    powerRTcorrStopSuccessFail2 = cell(1, numSessions);
    
    
    for iSession = 1 : numSessions

        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        numCh = length(sessionChannels);

        hilbert_sessionDir   = fullfile(subject_hilbertDir, sessionList{iSession});
        if ~exist(hilbert_sessionDir, 'dir'); continue; end
        
        powerRTcorr_sessionDir = fullfile(subject_powerRTcorrdir, sessionList{iSession});
        if ~exist(powerRTcorr_sessionDir, 'dir')
            mkdir(powerRTcorr_sessionDir);
        end
        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorrdir, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            mkdir(phaseRTcorr_sessionDir);
        end
        
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
        
        if ~exist(metadata_filename, 'file')
            error([metadata_filename ' could not be found.']);
        end
        load(metadata_filename);
        
        centerFreqs = mean(metadata.freqBands, 2);
        numFreqs = length(centerFreqs);
        
        powerRTcorrStopSuccessFail_metadata.freqs = centerFreqs;
        powerRTcorrStopSuccessFail_metadata.chNames = metadata.chNames;
        powerRTcorrStopSuccessFail_metadata.eventList = eventList;
        powerRTcorrStopSuccessFail_metadata.twin = twin;
        powerRTcorrStopSuccessFail_metadata.Fs = metadata.Fs;
        powerRTcorrStopSuccessFail_metadata.trialTypeSucc = trialType;
        powerRTcorrStopSuccessFail_metadata.trialTypeFail = trialtype2;
        powerRTcorrStopSuccessFail_metadata.completeChannels = 0;
        
        phaseRTcorrStopSuccessFail_metadata.freqs = centerFreqs;
        phaseRTcorrStopSuccessFail_metadata.chNames = metadata.chNames;
        phaseRTcorrStopSuccessFail_metadata.eventList = eventList;
        phaseRTcorrStopSuccessFail_metadata.twin = twin;
        phaseRTcorrStopSuccessFail_metadata.Fs = metadata.Fs;
        phaseRTcorrStopSuccessFail_metadata.trialTypeSucc = trialType;
        phaseRTcorrStopSuccessFail_metadata.trialTypeFail = trialtype2;
        phaseRTcorrStopSuccessFail_metadata.RTcutoffs = RTquantile_borders;
        
        numSamps = round(range(twin) * metadata.Fs);
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            hilbert_name = ['analytic_' sessionChannels{iCh}.name '.bin'];
            hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
            
            phase_RTcorr_name = ['phase_RTanalysisStopSuccessFail_' sessionChannels{iCh}.name '.mat'];
            phase_RTcorr_name = fullfile(phaseRTcorr_sessionDir, phase_RTcorr_name);
            
            power_RTcorr_name = ['power_RTcorrStopSuccessFail_' sessionChannels{iCh}.name '.mat'];
            power_RTcorr_name = fullfile(powerRTcorr_sessionDir, power_RTcorr_name);
            
            powerRTcorrStopSuccessFail = NaN(numEvents, numFreqs, numSamps);
             powerRTcorrStopSuccessFail2 = NaN(numEvents, numFreqs, numSamps);
            RTphasesStopSuccessFail    = cell(numEvents, numFreqs, numQuantiles);
            
            % load analytic signal around each event and calculate
            % correlation coefficients
            if exist(phase_RTcorr_name, 'file') && exist(power_RTcorr_name, 'file')
                continue
            end
            
            for iEvent = 1 : numEvents
                tic
                for iFreq = 1 : numFreqs
                    
                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    iFreq, ...
                                                    eventList{iEvent}, ...
                                                    trialType, ...
                                                    twin );
                    freqPower = abs(ansig) .^ 2;
                    freqPhase = angle(ansig);

                    ansig2 = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    iFreq, ...
                                                    eventList{iEvent}, ...
                                                    trialType2, ...
                                                    twin );
                    freqPower2 = abs(ansig2) .^ 2;
                    freqPhase2 = angle(ansig2);
                        
                        RTphasesStopSuccessFail{iEvent, iFreq, 1} = freqPhase(1, :);
                    RTphasesStopSuccessFail{iEvent, iFreq, 2} = freqPhase(1, :);

                    for iSamp = 1 : numSamps
                        [cc, p] = corr(freqPower(:, iSamp), RT{iSession}', ...
                                       'type','spearman');
                        powerRTcorrStopSuccessFail(iEvent, iFreq, iSamp) = cc;
                        
                        [cc2, p2] = corr(freqPower2(:, iSamp), RT2{iSession}', ...
                                       'type','spearman');
                        powerRTcorrStopSuccessFail2(iEvent, iFreq, iSamp) = cc2;
                    end
                        
                end
                toc
            end
            
            save(phase_RTcorr_name, 'RTphasesStopSuccessFail', 'phaseRTcorrStopSuccessFail_metadata');
            save(power_RTcorr_name, 'powerRTcorrStopSuccessFail', 'powerRTcorrStopSuccessFail_metadata','powerRTcorrStopSuccessFail2');
            
        end
    end
    
end
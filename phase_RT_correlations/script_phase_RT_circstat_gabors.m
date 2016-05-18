% script_phase_RT_circstat_gabors

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
gabor_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
% powerRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabors';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

% eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
% numEvents = length(eventList);
% twin = [-1 1];

trialType = 'correctgo';

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_gabor_dir = fullfile(gabor_directory, [implantID '_ps']);
%     subject_gabor_dir = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
%     subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    if ~exist(subject_gabor_dir, 'dir')
        disp([subject_gabor_dir ' not found. Skipping ' implantID '...'])
        continue
    end
%     if ~exist(subject_hilbertDir_025Hz, 'dir')
%         disp([subject_hilbertDir_025Hz ' not found. Skipping ' implantID '...'])
%         continue
%     end
    
    subject_phaseRTcorrdir = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr_gabors']);
    if ~exist(subject_phaseRTcorrdir, 'dir')
        mkdir(subject_phaseRTcorrdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    [RT, ~, sessionList] = collect_RT_MT_by_rat(channels, trialType);
    numSessions = length( sessionList );
    
    % establish RT quantiles for phase analysis
    allRT = RT{1};
    if numSessions > 1
        for iSession = 2 : numSessions
            allRT = [allRT, RT{iSession}];
        end
    end
%     allRT = sort(allRT);
%     numRT = length(allRT);
    
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

        gabor_sessionDir = fullfile(subject_gabor_dir, [sessionList{iSession} '_scalograms']);
%         hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
        if ~exist(gabor_sessionDir, 'dir'); continue; end
%         if ~exist(hilbert_sessionDir_025Hz, 'dir'); continue; end

        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorrdir, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            mkdir(phaseRTcorr_sessionDir);
        end
        
%         gabor_fn = [sessionList{iSession} '_' trialType '_scalograms.mat'];
%         gabor_fn = fullfile(gabor_sessionDir, gabor_fn);
%         if ~exist(gabor_fn, 'file')
%             error([gabor_fn ' could not be found.']);
%         end
%         metadata_filename_gabor = [sessionList{iSession} 'gabor_metadata.mat'];
%         metadata_filename_gabor = fullfile(gabor_sessionDir, metadata_filename_gabor);
%         metadata_filename_025Hz = [sessionList{iSession} 'hilbert_metadata.mat'];
%         metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename_025Hz);
%         if ~exist(metadata_filename_gabor, 'file')
%             error([metadata_filename_gabor ' could not be found.']);
%         end
%         if ~exist(metadata_filename_025Hz, 'file')
%             error([metadata_filename_025Hz ' could not be found.']);
%         end
%         md_1Hz   = load(metadata_filename_1Hz);
%         md_025Hz = load(metadata_filename_025Hz);
        
%         centerFreqs = mean(md_1Hz.metadata.freqBands, 2);
%         centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
% 
%         max025Hz = max(centerFreqs_025Hz);
%         freqList = [centerFreqs_025Hz; centerFreqs_1Hz(centerFreqs_1Hz > max025Hz)];
%             
%         numFreqs = length(freqList);
        
        
%         phaseRTcorr_metadata.eventList = eventList;
%         phaseRTcorr_metadata.twin = twin;
%         phaseRTcorr_metadata.Fs = md_1Hz.metadata.Fs;
        phaseRTcorr_metadata.trialType = trialType;
        
%         numSamps = round(range(twin) * md_1Hz.metadata.Fs);
            
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            disp(sprintf('working on session %s, %d of %d; channel %s, %d of %d', ...
                sessionList{iSession}, iSession, numSessions, ch.name, iCh, numCh))
            phaseRTcorr_metadata.chName = ch.name;
            
            gabor_name = [sessionChannels{iCh}.name '_' trialType '_scalograms.mat'];
            gabor_name = fullfile(gabor_sessionDir, gabor_name);
%             hilbert_name_025Hz = ['analytic_' sessionChannels{iCh}.name '.bin'];
%             hilbert_name_025Hz = fullfile(hilbert_sessionDir_025Hz, hilbert_name_025Hz);
            if ~exist(gabor_name,'file');continue;end
            
            phase_RTcorr_name = ['phase_RT_analysis_gabor_' sessionChannels{iCh}.name '.mat'];
            phase_RTcorr_name = fullfile(phaseRTcorr_sessionDir, phase_RTcorr_name);
            
            % load analytic signal around each event and calculate
            % correlation coefficients
            if exist(phase_RTcorr_name, 'file')
                continue
            end
            
            load(gabor_name);
            phaseRTcorr_metadata.freqs = scalogram_metadata.f;
            phaseRTcorr_metadata.twin = scalogram_metadata.twin;
            phaseRTcorr_metadata.t = scalogram_metadata.t;
            phaseRTcorr_metadata.eventList = scalogram_metadata.eventList;
            phaseRTcorr_metadata.Fs = scalogram_metadata.Fs;            
            
            numEvents = length(phaseRTcorr_metadata.eventList);
            numSamps = size(W,2);
            numFreqs = length(phaseRTcorr_metadata.freqs);

            circRTcorr = NaN(numEvents, numFreqs, numSamps);
            circRT_p = NaN(numEvents, numFreqs, numSamps);
        
            for iEvent = 1 : numEvents
                iEvent
                tic
                for iFreq = 1 : numFreqs
%                     iFreq
                    
%                     if freqList(iFreq) <= max025Hz
%                         hilbertDir = hilbert_025Hz_directory;
%                         freqIdx = find(centerFreqs_025Hz == freqList(iFreq));
%                     else
%                         hilbertDir = hilbert_1Hz_directory;
%                         freqIdx = find(centerFreqs_1Hz == freqList(iFreq));
%                     end

%                     ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
%                                                     freqIdx, ...
%                                                     eventList{iEvent}, ...
%                                                     trialType, ...
%                                                     twin, ...
%                                                     'hilbertdir', hilbertDir);
                    ansig = squeeze(W(iEvent,:,:,iFreq));
                    freqPhase = angle(ansig)';
                    
                    % GENERATE SURROGATE DISTRIBUTIONS? ALSO, ADD IN A
                    % CALCULATION WHERE ALL RTs ARE TESTED ACROSS SESSIONS?
                    
                    trialEventParams = getTrialEventParams('correctgo');
                    correctGOtrials = extractTrials2(ch.trials,trialEventParams);
                    % exclude trials that occur too early to
                    % get full gabor windows
                    if correctGOtrials(1).timestamps.cueOn < 2
                        current_RT = RT{iSession}(2:end);
                    else
                        current_RT = RT{iSession};
                    end
                    for iSamp = 1 : numSamps
%                         [circRTcorr(iEvent, iFreq, iSamp), ...
%                             circRT_p(iEvent, iFreq, iSamp)] = ...
%                             circ_corrcl(freqPhase(:, iSamp), RT{iSession}');
                        [circRTcorr(iEvent, iFreq, iSamp), ...
                            circRT_p(iEvent, iFreq, iSamp)] = ...
                            circ_corrcl(freqPhase(:, iSamp), current_RT');
                    end  
                end
                toc
            end
            
            save(phase_RTcorr_name, 'circRTcorr', 'circRT_p', 'phaseRTcorr_metadata');
            
        end
    end
    
end
% script_analyze_STOPsuccvfail_phases

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

nBoot = 1000;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        disp([subject_stopPhaseDir ' not found. Skipping ' implantID '...'])
        continue
    end

    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.task = 3;    % only stop-signal sessions
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numCh = length(sessionChannels);
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end
        
        for iCh = 1 : numCh
            iCh
            
            STOPanal_metadata.eventsComplete = 0;
            STOPanal_metadata_loaded = false;
            
            ch = sessionChannels{iCh};
            wwtest_saveName = ['STOPsucc_fail_wwtest_' ch.name '.mat'];
            wwtest_saveName = fullfile(session_stopPhaseDir, wwtest_saveName);
            bootstrap_saveName = ['STOPsucc_fail_boot_' ch.name '.mat'];
            bootstrap_saveName = fullfile(session_stopPhaseDir, bootstrap_saveName);
            if exist(wwtest_saveName, 'file'); 
                load(wwtest_saveName)
                if (STOPanal_metadata.eventsComplete == length(STOPanal_metadata.eventList)) && ...
                   (STOPanal_metadata.freqComplete == size(STOPanal_metadata.freqBands, 1))
                    continue;
                end
                STOPanal_metadata_loaded = true;
                load(bootstrap_saveName);
            end
            
            stopPhase_name = ['stopPhases_' ch.name '.mat'];
            stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
            if ~exist(stopPhase_name, 'file'); continue; end
            
            load(stopPhase_name);
            STOPanal_metadata.Fs = STOPmetadata.Fs;
            STOPanal_metadata.freqBands = STOPmetadata.freqBands;
            STOPanal_metadata.eventList = STOPmetadata.eventList;
            STOPanal_metadata.chName = STOPmetadata.chName;
            STOPanal_metadata.twin = STOPmetadata.twin;
            
            numEvents = length(STOPmetadata.eventList);
            numFreqs  = size(STOPmetadata.freqBands, 1);
            numSamps  = size(correctSTOP_phase, 4);
            
            if ~STOPanal_metadata_loaded
                ww_p   = zeros(numEvents, numFreqs, numSamps);       % wwtest p-value
                boot_p = zeros(numEvents, numFreqs, numSamps);       % bootstrapped p-value
                r      = zeros(2, numEvents, numFreqs, numSamps);    % mean resultant vector length
                theta  = zeros(2, numEvents, numFreqs, numSamps);    % mean resultant vector phase
            end
            
            for iEvent = STOPanal_metadata.eventsComplete + 1 : numEvents
                tic
                if ~STOPanal_metadata_loaded
                    STOPanal_metadata.freqComplete = 0;
                end
                for iFreq = STOPanal_metadata.freqComplete + 1 : numFreqs
                    iFreq
                    for i_t = 1 : numSamps
                        
                        phases_correct = squeeze(correctSTOP_phase(iEvent, iFreq, :, i_t));
                        phases_failed  = squeeze(failedSTOP_phase(iEvent, iFreq, :, i_t));
                        [ww_p(iEvent, iFreq, i_t), ~] = circ_wwtest(phases_correct, phases_failed);
                        temp = boot_angle_test(phases_correct, phases_failed, 'nboot', nBoot);
                        boot_p(iEvent, iFreq, i_t) = boot_angle_test(phases_correct, phases_failed, 'nboot', nBoot);
                        
                        r(1, iEvent, iFreq, i_t) = abs(mean(exp(1i * phases_correct)));
                        theta(1, iEvent, iFreq, i_t) = angle(mean(exp(1i * phases_correct)));
                        
                        r(2, iEvent, iFreq, i_t) = abs(mean(exp(1i * phases_failed)));
                        theta(2, iEvent, iFreq, i_t) = angle(mean(exp(1i * phases_failed)));
                        
                        
                    end
                    STOPanal_metadata.freqComplete = STOPanal_metadata.freqComplete + 1;
                    save(wwtest_saveName, 'ww_p', 'r', 'theta', 'STOPanal_metadata');
                    save(bootstrap_saveName, 'boot_p', 'r', 'theta', 'STOPanal_metadata');
                end
                STOPanal_metadata_loaded = false;
                toc
            end    % for iEvent...
            
        end
        
    end
    
end
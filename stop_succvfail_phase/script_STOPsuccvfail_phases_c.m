% script_STOPsuccvfail_phases

% script to go through all STOP sessions and save continuous phase values
% for stop-success and stop-fail trials

% chDB_directory    = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
% phase_stop_succvfail_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
% power_stop_succvfail_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
numEventTypes = length(eventList);
twin = [-1 1];

for i_chDB = 3 : 3%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
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
    
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        mkdir(subject_stopPhaseDir);
    end
    subject_stopPowerDir = fullfile(power_stop_succvfail_directory, [implantID '_stopPower']);
    if ~exist(subject_stopPowerDir, 'dir')
        mkdir(subject_stopPowerDir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    cp = initChanParams();
    cp.task = 3;
    chList = extractChannels(cp, channels);
    channels = channels(chList);
    
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
            mkdir(session_stopPhaseDir);
        end
        session_stopPowerDir = fullfile(subject_stopPowerDir, sessionList{iSession});
        if ~exist(session_stopPowerDir, 'dir')
            mkdir(session_stopPowerDir);
        end
        hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
        hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
        if ~exist(hilbert_sessionDir_1Hz, 'dir')
            disp([hilbert_sessionDir_1Hz ' could not be found.']);
            continue;
        end
        if ~exist(hilbert_sessionDir_025Hz, 'dir')
            disp([hilbert_sessionDir_025Hz ' could not be found.']);
            continue;
        end
        
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename_1Hz = fullfile(hilbert_sessionDir_1Hz, metadata_filename);
        metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename);

        if ~exist(metadata_filename_1Hz, 'file')
            continue
%             error([metadata_filename ' could not be found.']);
        end
        if ~exist(metadata_filename_025Hz, 'file')
            continue
%             error([metadata_filename ' could not be found.']);
        end
        md_1Hz   = load(metadata_filename_1Hz);
        md_025Hz = load(metadata_filename_025Hz);
        
        centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
        centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
        max_025Hz_freq = max(centerFreqs_025Hz);
        centerFreqs = [centerFreqs_025Hz; ...
                       centerFreqs_1Hz];
        numFreqs = length(centerFreqs);
        
        numSamps = round(range(twin) * md_1Hz.metadata.Fs);

        STOPmetadata.Fs = md_1Hz.metadata.Fs;
        STOPmetadata.freqBands = [md_025Hz.metadata.freqBands; ...
                                  md_1Hz.metadata.freqBands];
        STOPmetadata.eventList = eventList;
        STOPmetadata.twin      = twin;
            
        ch = sessionChannels{1};
        if ch.task ~= 3; continue; end    % only stop sessions
        
        trialEventParams = getTrialEventParams('correctstop');
        trialList_correct = extractTrials(ch.trials, trialEventParams);
        num_correctSTOP = length(trialList_correct);
        
        trialEventParams = getTrialEventParams('failedstop');
        trialList_failed = extractTrials(ch.trials, trialEventParams);
        num_failedSTOP = length(trialList_failed);
        
        for iCh = 1 : numCh
            iCh
            
            ch = sessionChannels{iCh};
            saveName_phase = ['stopPhases_' ch.name '.mat'];
            saveName_phase = fullfile(session_stopPhaseDir, saveName_phase);
            if exist(saveName_phase, 'file'); continue; end
            
            saveName_power = ['stopPower_' ch.name '.mat'];
            saveName_power = fullfile(session_stopPowerDir, saveName_power);
            
%             metadata_filename = [ch.session 'hilbert_metadata.mat'];
%             metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
%             if ~exist(metadata_filename, 'file')
%                 error('metadata file not found.');
%                 continue;
%             end


            STOPmetadata.chName = ch.name;
            
            failedSTOP_power = zeros(numEventTypes, numFreqs, num_failedSTOP, numSamps);
            correctSTOP_power = zeros(numEventTypes, numFreqs, num_correctSTOP, numSamps);

            failedSTOP_phase = zeros(numEventTypes, numFreqs, num_failedSTOP, numSamps);
            correctSTOP_phase = zeros(numEventTypes, numFreqs, num_correctSTOP, numSamps);
            
%             hilbert_name = ['analytic_' ch.name '.bin'];
%             hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);

            for iEvent = 1 : numEventTypes
                iEvent
                tic
                for iFreq = 1 : numFreqs

                    if centerFreqs(iFreq) <= max_025Hz_freq
                        hilbert_directory = hilbert_025Hz_directory;
                        freq_idx = find(centerFreqs_025Hz == centerFreqs(iFreq));
                    else
                        hilbert_directory = hilbert_1Hz_directory;
                        freq_idx = find(centerFreqs_1Hz == centerFreqs(iFreq));
                    end
                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    freq_idx, ...
                                                    eventList{iEvent}, ...
                                                    'correctstop', ...
                                                    twin, ...
                                                    'hilbertdir', hilbert_directory );
                                                
                    correctSTOP_power(iEvent, iFreq, :, :) = abs(ansig) .^ 2;
                    correctSTOP_phase(iEvent, iFreq, :, :) = angle(ansig);

                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    freq_idx, ...
                                                    eventList{iEvent}, ...
                                                    'failedstop', ...
                                                    twin, ...
                                                    'hilbertdir', hilbert_directory );
                                                
                    failedSTOP_power(iEvent, iFreq, :, :) = abs(ansig) .^ 2;
                    failedSTOP_phase(iEvent, iFreq, :, :) = angle(ansig);
                    
                end    % for iFreq...
                toc
            end    % for iEvent
            
            save(saveName_phase, 'correctSTOP_phase', 'failedSTOP_phase', 'STOPmetadata');
            save(saveName_power, 'correctSTOP_power', 'failedSTOP_power', 'STOPmetadata');
            
        end    % for iCh...
        
    end    % for iSession...
    
end    % for i_chDB
                    

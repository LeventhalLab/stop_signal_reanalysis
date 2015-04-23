% script to update the power calculations for correctgo and wronggo trials
% to include the nose out event.

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

powerSpectrogramDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectrograms';
% phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventList{1} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideOut'};
eventList{2} = eventList{1};

trialTypeList = {'correctgo', 'wronggo'};

numTrialTypes = length(trialTypeList);

eventtWin = zeros(numTrialTypes, 2);
for iTrialType = 1:numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
end

for i_chDB = 2 : length(chDB_list)
    
    if i_chDB > 7
        hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
    else
        hilbert_025Hz_directory = '/Volumes/RecordingsLeventhal2/stop-sig_reanalysis BU/Hilbert transformed LFP 025 Hz bins';
    end

    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
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
    
    subject_powerSpectDir = fullfile(powerSpectrogramDir, [implantID '_powerSpectrograms']);
    if ~exist(subject_powerSpectDir, 'dir')
        mkdir(subject_powerSpectDir);
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    if i_chDB == 2 && iTrialType == 1
        startSession = 21;
    else
        startSession = 21;
    end
    for iTrialType = 1 : length(trialTypeList)
        
        trialType = trialTypeList{iTrialType}
        numEvents = length(eventList{iTrialType});
        
        for iSession = startSession : numSessions
            
            powerSpect_sessionDir = fullfile(subject_powerSpectDir, sessionList{iSession});
            if ~exist(powerSpect_sessionDir, 'dir')
                mkdir(powerSpect_sessionDir);
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.task = -1;
            
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
            
            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            
            cp = initChanParams();
            cp.tetrode = {'e2', 'e3', 'e03','e03'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels);continue;end
        
            numCh = length(sessionChannels);
            
            hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
            hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
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
            
            max025Hz = max(centerFreqs_025Hz);
            num_freq_025 = length(centerFreqs_025Hz);
            
            ch = sessionChannels{1};    % just to get the trial list correct
            trialEventParams = getTrialEventParams(trialType);
            trList = extractTrials(ch.trials, trialEventParams);

            % figure out how many VALID trials there are (that is, trials that don't
            % occur to close to the beginning or end of the session so that we can't
            % extract a full window around all events
            numTrials = length(trList);
            minValidTrials = numTrials;
            for iEventType = 1 : numEvents
                eventName = eventList{iTrialType}{iEventType};
                numValidTrials = 1;
                for iTr = 1 : numTrials
                    if ~isfield(ch.trials(trList(numValidTrials)).timestamps, eventName)
                        if numValidTrials == 1
                            trList = trList(2:end,:);
                        elseif numValidTrials == length(trList)
                            trList = trList(1:numValidTrials);
                        else
                            trList = [trList(1:numValidTrials-1); trList(numValidTrials+1:end)];
                        end
                        continue
                    end
                    event_ts = ch.trials(trList(numValidTrials)).timestamps.(eventName);

                    tlim = event_ts + eventtWin(iTrialType,:);
                    % make sure the time window doesn't extend before the start of the
                    % recording or after the end of the recording
                    if any(tlim < 0) || any(tlim > md_1Hz.metadata.duration)
                        if numValidTrials == 1
                            trList = trList(2:end,:);
                        elseif numValidTrials == length(trList)
                            trList = trList(1:numValidTrials);
                        else
                            trList = [trList(1:numValidTrials-1); trList(numValidTrials+1:end)];
                        end
                        continue;
                    end
                    numValidTrials = numValidTrials + 1;
                end
                numTrials = numValidTrials - 1;
                if numValidTrials - 1 < minValidTrials
                    minValidTrials = numValidTrials - 1;
                end
            end
            numValidTrials = minValidTrials;
            
            
            
            
            for iCh = 1 : numCh
                disp(sprintf('trialtype: %s, session %s, %d of %d; channel %d of %d', ...
                    trialType, sessionList{iSession}, iSession, numSessions, iCh, numCh))

                ch = sessionChannels{iCh};
            
                powerSpect_name = ['powerSpect_' trialType '_' sessionChannels{iCh}.name '.mat'];
                powerSpect_name = fullfile(powerSpect_sessionDir, powerSpect_name);
                
                if exist(powerSpect_name, 'file')
                    load(powerSpect_name)
                else
                    continue;
                end
                
                if strcmpi(powerSpect_metadata.eventList{end}, eventList{iTrialType}{end})
                    continue
                end
                % 'noseSideOut' isn't the last event, so we have some new
                % calculations to do
                
                new_powerSpect = zeros(length(eventList{iTrialType}), size(powerSpect,2),size(powerSpect,3));
                new_powerSpect(1:length(eventList{iTrialType})-1,:,:) = powerSpect;
                
                powerSpect_metadata.eventList = eventList{iTrialType};
                f = powerSpect_metadata.freqs;
                numFreqs = length(f);
                
                twin = powerSpect_metadata.eventtWin;
                
                for iEvent = length(eventList{iTrialType}):length(eventList{iTrialType})
                    
                    for iFreq = 1 : numFreqs
                        
                        if iFreq <= num_freq_025
                            activeHilbertDir = hilbert_025Hz_directory;
                            freqIdx = iFreq;
                        else
                            activeHilbertDir = hilbert_1Hz_directory;
                            freqIdx = find(centerFreqs_1Hz == f(iFreq));
                        end
                        
                        ansig = getAnalyticAroundEvent_20140916( ch, ...
                                                                 trList, ...
                                                                 freqIdx, ...
                                                                 eventList{iTrialType}{iEvent}, ...
                                                                 twin, ...
                                                                 'hilbertdir', activeHilbertDir);
                                                             
                        freqPower = abs(ansig) .^ 2;
                        new_powerSpect(iEvent, iFreq, :) = squeeze(mean(freqPower,1));
                        
                    end    % for iFreq...
                    
                end    % for iEvent...
                
                powerSpect = new_powerSpect;
                save(powerSpect_name, 'powerSpect', 'powerSpect_metadata');
                
            end    % for iCh...
            
        end    % for iSession...
        
    end    % for iTrialType...
    
end    % for i_chDB...
                
                
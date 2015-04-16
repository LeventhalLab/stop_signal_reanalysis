% script_phaseAmp_surrogates_byTrial_t_mrl_20140324c3

% rev history
%  3/24/2013 - added ability to include phase frequencies from the 1 Hz
%       spacing analytic signal files.
bitOrder = 'b';
low_freq_range  = [0 21];
high_freq_range = [10 101];

% in case we need to go back to older versions and add on
% old_low_freqs = [1.5000
%     1.7500
%     2.0000
%     2.2500
%     2.5000
%     2.7500
%     3.0000
%     3.2500
%     3.5000
%     3.7500
%     4.0000
%     4.2500
%     4.5000
%     4.7500
%     5.0000
%     5.2500
%     5.5000
%     5.7500
%     6.0000
%     6.2500
%     6.5000
%     6.7500
%     7.0000
%     7.2500
%     7.5000
%     7.7500
%     8.0000
%     8.2500
%     8.5000
%     8.7500
%     9.0000
%     9.2500
%     9.5000
%     9.7500
%    10.0000];

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventLists{1} = {'noseCenterIn'};
eventLists{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
eventLists{3} = eventLists{2};
eventLists{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventLists{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventLists{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventLists{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

numTrialTypes = length(trialTypeList);
eventtWin   = zeros(numTrialTypes, 2);
eventtWin(1,:) = [-1 2];   % for analysis of all trials
for iTrialType = 2 : numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
end

numSurrogates = 200;
maxSkip = 3;    % in seconds
minSkip = 0;

for i_chDB = 4 : 4%length(chDB_list)
    
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

    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );

    if i_chDB==4
        startTrialType = 2;
    else
        startTrialType = 2;
    end
    for iTrialType = startTrialType : 2%length(trialTypeList)
        trialType = trialTypeList{iTrialType};
        eventList = eventLists{iTrialType};
        numEventTypes = length(eventList);
        
        surrogate_phaseAmp_metadata.low_freq_range  = low_freq_range;
        surrogate_phaseAmp_metadata.high_freq_range = high_freq_range;
        surrogate_phaseAmp_metadata.eventList       = eventList;
        surrogate_phaseAmp_metadata.trialType       = trialType;
        surrogate_phaseAmp_metadata.eventtWin       = eventtWin(iTrialType, :);
%         surrogate_phaseAmp_metadata.analysisWin     = analysisWin(iTrialType);
%         surrogate_phaseAmp_metadata.stepSize        = stepSize(iTrialType);
        
        for iSession = 11 : 11
            
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                mkdir(phaseAmp_sessionDir);
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            
            if ~isempty(strfind(trialType, 'nogo'))
                cp.task = 4;
            elseif ~isempty(strfind(trialType, 'stop'))
                cp.task = 3;
            else
                cp.task = -1;
            end
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
            if isempty(sessionChannels);continue;end

            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels); continue; end
            
            ch = sessionChannels{1};
            
            trialEventParams = getTrialEventParams(trialType);
            trList = extractTrials(ch.trials, trialEventParams);
            numTrials = length(trList);
            

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
            
            % figure out how many VALID trials there are (that is, trials that don't
            % occur to close to the beginning or end of the session so that we can't
            % extract a full window around all events
            
            minValidTrials = numTrials;
            for iEventType = 1 : numEventTypes
                eventName = eventList{iEventType};
                numValidTrials = 1;
                for iTr = 1 : numTrials
                    if ~isfield(ch.trials(trList(numValidTrials)).timestamps, eventName); continue; end
                    event_ts = ch.trials(trList(numValidTrials)).timestamps.(eventName);

                    tlim = event_ts + surrogate_phaseAmp_metadata.eventtWin;
                    % make sure the time window doesn't extend before the start of the
                    % recording or after the end of the recording
                    if any(tlim < 0) || any(tlim > md_1Hz.metadata.duration)
                        if numValidTrials == 1
                            trList = trList(2:end,:);
                        elseif numValidTrials == length(trList)
                            trList = trList(1:numValidTrials);
                        else
                            trList = [trList(1:numValidTrials-1); trList(numValidTrials:end)];
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
            
            centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
            centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
            
            max025Hz = max(centerFreqs_025Hz);

            low_freq_idx_025  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
            low_freq_idx_1    = find(centerFreqs_1Hz >= max025Hz & centerFreqs_1Hz <= low_freq_range(2));
            low_freq_idx      = [low_freq_idx_025; low_freq_idx_1];
            high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

            low_freqs  = [centerFreqs_025Hz(low_freq_idx_025); centerFreqs_1Hz(low_freq_idx_1)];
            high_freqs = centerFreqs_1Hz(high_freq_idx);
            
            surrogate_phaseAmp_metadata.low_freq = low_freqs;
            surrogate_phaseAmp_metadata.high_freq = high_freqs;
            
            num_low_freq_025 = length(low_freq_idx_025);
            num_low_freq  = length(low_freq_idx);
            num_high_freq = length(high_freq_idx);

            surrogate_phaseAmp_metadata.Fs = md_025Hz.metadata.Fs;
%             surrogate_phaseAmp_metadata.chList = cell(1, numCh);
%             surrogate_phaseAmp_metadata.regionList = cell(1, numCh);
            
            sampsPerTrial = round(range(eventtWin(iTrialType, :)) * surrogate_phaseAmp_metadata.Fs);
            totSamps = sampsPerTrial * numTrials;
            
            minSkipSamps = round(minSkip * surrogate_phaseAmp_metadata.Fs);
            maxSkipSamps = round(maxSkip * surrogate_phaseAmp_metadata.Fs);
            skipSampRange = maxSkipSamps - minSkipSamps;
            skip = minSkipSamps + ceil(skipSampRange .* rand(numSurrogates,1));
%             skip(find(skip>maxSkip))=[];
%             skip(find(skip<minSkip))=[];
            
            for iCh = 1: numCh
                ch = sessionChannels{iCh};
                
                disp(sprintf('%s, %s, %d of %d sessions, %d of %d channels in session %s', ...
                    trialType, ...
                    ch.name, ...
                    iSession, ...
                    numSessions, ...
                    iCh, ...
                    numCh, ...
                    sessionList{iSession}));
                
                phaseAmp_name_mrl = [ch.name '_surrTrial_t_' trialType '_phase_amp_mrl.mat'];
                phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
                if exist(phaseAmp_name_mrl, 'file')
                    load(phaseAmp_name_mrl);
                    if ~isfield(surrogate_phaseAmp_metadata, 'numEventsComplete') && ...
                        (iTrialType == 1)
                        if surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq
                            surrogate_phaseAmp_metadata.numEventsComplete = 1;
                        else
                            surrogate_phaseAmp_metadata.numEventsComplete = 0;
                        end
                    end
                    
                    if (surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq) && ...
                       (surrogate_phaseAmp_metadata.numEventsComplete == numEventTypes)
                        continue
                    end
                    numLowFreqsComplete = surrogate_phaseAmp_metadata.numLowFreqsComplete;
                    if isfield(surrogate_phaseAmp_metadata, 'numEventsComplete')
                        numEventsComplete = surrogate_phaseAmp_metadata.numEventsComplete;
                    end
                    surrogate_phaseAmp_metadata.low_freq = low_freqs;
                else
                    numLowFreqsComplete = 0;
                    numEventsComplete = 0;
                    surrogate_mean = zeros(numEventTypes, num_low_freq, num_high_freq, sampsPerTrial);
                    surrogate_std  = zeros(numEventTypes, num_low_freq, num_high_freq, sampsPerTrial);
                end
                
                surrogate_phaseAmp_metadata.chName = ch.name;
                surrogate_phaseAmp_metadata.region = ch.location.name;
                
                surrogate_phaseAmp_metadata.numEventsComplete = numEventsComplete;
                for iEventType = numEventsComplete + 1 : numEventTypes
                    iEventType
                    
                    as2 = zeros(num_high_freq, numValidTrials, sampsPerTrial);
                    as1 = zeros(num_low_freq, numValidTrials, sampsPerTrial);

                    for i_f2 = 1 : num_high_freq

                        temp = getAnalyticAroundEvent_20140916( ch, ...
                                                                trList, ...
                                                                high_freq_idx(i_f2), ...
                                                                eventList{iEventType}, ...
                                                                eventtWin(iTrialType, :), ...
                                                                'hilbertdir', hilbert_1Hz_directory );
                        as2(i_f2, :, :) = temp;
                    end
                    
%                     as2 = squeeze(reshape(as2, num_high_freq, numTrials * sampsPerTrial, 1));
                    sig_amp = abs(as2);

                    
                    for i_f1 = numLowFreqsComplete + 1 : num_low_freq
                        
                        if i_f1 <= num_low_freq_025 
                            activeHilbertDir = hilbert_025Hz_directory;
                        else
                            activeHilbertDir = hilbert_1Hz_directory;
                        end

                        temp = getAnalyticAroundEvent_20140916( ch, ...
                                                                trList, ...
                                                                low_freq_idx(i_f1), ...
                                                                eventList{iEventType}, ...
                                                                eventtWin(iTrialType, :), ...
                                                                'hilbertdir', activeHilbertDir );
                                                            
                        as1(i_f1, :, :) = temp;

                    end
                    
%                     as1 = squeeze(reshape(as1, num_high_freq, numTrials * sampsPerTrial, 1));
                    phase_angles = angle(as1);
                    for i_f1 = numLowFreqsComplete + 1 : num_low_freq
                        i_f1
                        for i_f2 = 1 : num_high_freq

                            if low_freqs(i_f1) >= high_freqs(i_f2); continue; end

                            amp_f2 = squeeze(sig_amp(i_f2, :, :));
                            surrogate_mrl = zeros(numSurrogates, sampsPerTrial);
                            for iSurrogate = 1 : numSurrogates
%                                 surrogate_amp(iSurrogate, :) = [squeeze(sig_amp(i_f2, :, skip(iSurrogate):end)), ...
%                                                                    squeeze(sig_amp(i_f2, :, 1:skip(iSurrogate)-1))];
                                surrogate_amp = circshift(amp_f2, [0, skip(iSurrogate)]);
                                surrogate_phaseAmp = surrogate_amp.*exp(1i*squeeze(phase_angles(i_f1, :, :)));
                                meanTrial = mean(surrogate_phaseAmp, 1);
                                surrogate_mrl(iSurrogate, :) = abs(meanTrial);              % SHOULD THIS BE ABS?
%                                 surrogate_mrl(iSurrogate) = abs(mean(surrogate_amp.*exp(1i*squeeze(phase_angles(i_f1, :)))));
                            end    % for iSurrogate...
                            surrogate_mean(iEventType, i_f1, i_f2, :) = mean(surrogate_mrl, 1);
                            surrogate_std(iEventType, i_f1, i_f2, :) = std(surrogate_mrl, 0 , 1);

                        end    % for i_f2...
                        
                        numLowFreqsComplete = numLowFreqsComplete + 1;
                        surrogate_phaseAmp_metadata.numLowFreqsComplete = numLowFreqsComplete;
                        save(phaseAmp_name_mrl, 'surrogate_mean', 'surrogate_std', 'surrogate_phaseAmp_metadata');

                    end    % for i_f1...
                
                    numEventsComplete = numEventsComplete + 1;
                    surrogate_phaseAmp_metadata.numEventsComplete = numEventsComplete;
                    numLowFreqsComplete = 0;
                    
                end    % for iEventType
                
                save(phaseAmp_name_mrl, 'surrogate_mean', 'surrogate_std', 'surrogate_phaseAmp_metadata');
                
                
            end    % for iCh = 1 : numCh
        end    % for iSession...
        
    end    % for iTrialType...
end    % for i_chDB...
                
% script to update the phase amplitude coupling surrogate distribution 
% calculations for correctgo and wronggo trials to include the nose side out event.


% AFTER THIS SCRIPT IS RUN ON A DATA SET, AND script_update_phaseAmp_calcs
% IS RUN ON THE SAME SET, NEED TO RERUN script_plot_windowed_phaseAmp_z_byRegion_20140404
% FOR THAT SET, FOLLOWED BY script_review_windowed_phaseAmp_acrossSessions_20140430

bitOrder = 'b';
low_freq_range  = [0 21];
high_freq_range = [10 101];

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/RecordingsLeventhal2/stop-sig_reanalysis BU/Hilbert transformed LFP 025 Hz bins';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventLists{1} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideOut'};
eventLists{2} = eventLists{1};

trialTypeList = {'correctgo', 'wronggo'};

numTrialTypes = length(trialTypeList);

eventtWin = zeros(numTrialTypes, 2);
for iTrialType = 1:numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
end

numSurrogates = 200;
maxSkip = 3;    % in seconds
minSkip = 0;

for i_chDB = 1 : length(chDB_list)
    
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
    
    for iTrialType = 1 : length(trialTypeList)
        
        trialType = trialTypeList{iTrialType};
        eventList = eventLists{iTrialType};
        numEventTypes = length(eventList);
        
        for iSession = 1:numSessions
            
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                continue;
            end
            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.task = -1;

            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
            if isempty(sessionChannels);continue;end

            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels); continue; end
            
            trialEventParams = getTrialEventParams(trialType);
            trList = extractTrials(sessionChannels{1}.trials, trialEventParams);
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
            
            for iCh = 1 : numCh
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
                    if ~isfield(surrogate_phaseAmp_metadata, 'numEventsComplete')
                        if iTrialType == 1
                            if surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq
                                surrogate_phaseAmp_metadata.numEventsComplete = 1;
                            else
                                surrogate_phaseAmp_metadata.numEventsComplete = 0;
                            end
                        else
                            if surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq
                                surrogate_phaseAmp_metadata.numEventsComplete = 1;
                                surrogate_phaseAmp_metadata.numLowFreqsComplete = 0;
                            else
                                surrogate_phaseAmp_metadata.numEventsComplete = 0;
                            end
                        end        
                    end
                    
                    if (surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq) && ...
                       (surrogate_phaseAmp_metadata.numEventsComplete == numEventTypes) && ...
                       (strcmpi(surrogate_phaseAmp_metadata.eventList{end},eventList{end}))
                        continue
                    end
                    numLowFreqsComplete = surrogate_phaseAmp_metadata.numLowFreqsComplete;
                    if isfield(surrogate_phaseAmp_metadata, 'numEventsComplete')
                        numEventsComplete = surrogate_phaseAmp_metadata.numEventsComplete;
                    end
                    surrogate_phaseAmp_metadata.low_freq = low_freqs;
                    
                    if ~strcmpi(surrogate_phaseAmp_metadata.eventList{end},eventList{end})
                        surrogate_phaseAmp_metadata.numLowFreqsComplete = 0;
                    end
                else
                    continue;
                end
                    
                if size(surrogate_mean, 1) < length(eventList)
                    new_surrogate_mean = zeros(length(eventList),...
                                               size(surrogate_mean,2), ...
                                               size(surrogate_mean,3), ...
                                               size(surrogate_mean,4));

                    new_surrogate_std  = zeros(length(eventList),...
                                               size(surrogate_mean,2), ...
                                               size(surrogate_mean,3), ...
                                               size(surrogate_mean,4));

                    new_surrogate_mean(1:length(eventList)-1,:,:,:) = surrogate_mean;
                    new_surrogate_std(1:length(eventList)-1,:,:,:) = surrogate_std;

                    surrogate_mean = new_surrogate_mean;
                    surrogate_std = new_surrogate_std;
                end

                for iEventType = length(eventList) : length(eventList)
                    iEventType
                    
                    as2 = zeros(num_high_freq, numTrials, sampsPerTrial);
                    as1 = zeros(num_low_freq, numTrials, sampsPerTrial);

                    for i_f2 = 1 : num_high_freq

                        temp = getAnalyticAroundEvent_20140221( ch, ...
                                                               high_freq_idx(i_f2), ...
                                                               eventList{iEventType}, ...
                                                               trialType, ...
                                                               eventtWin(iTrialType, :), ...
                                                               'hilbertdir', hilbert_1Hz_directory );
                        as2(i_f2, :, :) = temp;
                    end
                    
                    sig_amp = abs(as2);
                    
                    for i_f1 = numLowFreqsComplete + 1 : num_low_freq
                        
                        if i_f1 <= num_low_freq_025 
                            activeHilbertDir = hilbert_025Hz_directory;
                        else
                            activeHilbertDir = hilbert_1Hz_directory;
                        end

                        temp = getAnalyticAroundEvent_20140221( ch, ...
                                                                low_freq_idx(i_f1), ...
                                                                eventList{iEventType}, ...
                                                                trialType, ...
                                                                eventtWin(iTrialType, :), ...
                                                                'hilbertdir', activeHilbertDir );
                                                            
                        as1(i_f1, :, :) = temp;

                    end
                    
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
                        surrogate_phaseAmp_metadata.eventList = eventList;
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
                                               
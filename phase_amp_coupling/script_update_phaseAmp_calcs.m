% script to update the phase amplitude coupling calculations for correctgo
% and wronggo trials to include the nose side out event.

% script to update the phase amplitude coupling surrogate distribution 
% calculations for correctgo and wronggo trials to include the nose side out event.


% test this Wednesday ******************************************



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

for i_chDB = 1 : length(chDB_list)
    
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

    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        continue;
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iTrialType = 1 : length(trialTypeList)
        trialType = trialTypeList{iTrialType}
        eventList = eventLists{iTrialType};
        
        if i_chDB == 1 && iTrialType == 1
            startSession = 4;
        else
            startSession = 1;
        end
        for iSession = startSession : numSessions
            
            disp(sprintf('trialtype: %s, session %s, %d of %d', trialType, sessionList{iSession}, iSession, numSessions))
            
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                continue;
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.task = -1;
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
            
            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            
            cp = initChanParams();
            cp.tetrode = {'e2', 'e3', 'e02','e03'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels); continue; end

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
            
            num_low_freq_025 = length(low_freq_idx_025);
            num_low_freq  = length(low_freq_idx);
            num_high_freq = length(high_freq_idx);
            
            for iCh = 1 : numCh             
                
                ch = sessionChannels{iCh};
                if any(ch.wire.markedGood) == 0; continue; end
                
                disp(sprintf('%s, %d of %d sessions, %d of %d channels in session %s', ...
                    ch.name, ...
                    iSession, ...
                    numSessions, ...
                    iCh, ...
                    numCh, ...
                    sessionList{iSession}));
                
                phaseAmp_name_mrl = [ch.name '_' trialType '_phase_amp_mrl.mat'];
                phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
                if ~exist(phaseAmp_name_mrl, 'file'); continue; end
                
                load(phaseAmp_name_mrl)
                
                if strcmpi(phaseAmp_metadata.eventList{end},eventList{end})
                    continue;
                end
                phaseAmp_metadata.eventList = eventList;
                new_re_mrv = re_mrv;
                new_im_mrv = im_mrv;
                re_mrv = zeros(length(eventList), size(re_mrv,2),size(re_mrv,3),size(re_mrv,4));
                im_mrv = zeros(length(eventList), size(re_mrv,2),size(re_mrv,3),size(re_mrv,4));
                re_mrv(1:length(eventList)-1,:,:,:) = new_re_mrv;
                im_mrv(1:length(eventList)-1,:,:,:) = new_im_mrv;
                
                [mrv, low_freqs, high_freqs] = phase_amp_Canolty_mrl_20140324(ch, ...
                                                                              eventList{end}, ...
                                                                              eventtWin(iTrialType, :), ...
                                                                              trialType, ...
                                                                              'lowfreqrange', low_freq_range, ...
                                                                              'highfreqrange', high_freq_range, ...
                                                                              'hilbert1Hzdir', hilbert_1Hz_directory, ...
                                                                              'hilbert025Hzdir', hilbert_025Hz_directory );
                re_mrv(length(eventList),:,:,:) = real(mrv);
                im_mrv(length(eventList),:,:,:) = imag(mrv);
                
                phaseAmp_metadata.low_freqs = low_freqs; phaseAmp_metadata.high_freqs = high_freqs;
                save(phaseAmp_name_mrl, 're_mrv', 'im_mrv', 'phaseAmp_metadata');
                
            end    % for iCh...
            
        end    % for iSession...
        
    end    % for iTrialType...
    
end    % for i_chDB...
                
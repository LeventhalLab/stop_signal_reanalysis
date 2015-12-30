% script to calculate scalograms on a per-trial basis only for tetrodes in
% areas we care about

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
scalogramDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
% lfp_root          = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs';
sample_lfp_name = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/high_cutoff_stop-singal LFPs/IM164_HF_LFPs/IM164_20091112_13-56-32.hsdf';
lfp_root          = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/high_cutoff_stop-singal LFPs';

twin_buffer = 1;    % add an extra second before and after each time window to eliminate edge effects

machineFormat = 'b';

header = getHSDHeader( sample_lfp_name );
Fs = lfpFs( header );
scalogram_metadata.Fs = Fs;
startSample = round(Fs * twin_buffer);

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

eventLists{1} = {'noseCenterIn'};
eventLists{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideOut'};
eventLists{3} = eventLists{2};
eventLists{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventLists{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventLists{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventLists{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

numTrialTypes = length(trialTypeList);
eventtWin   = zeros(numTrialTypes, 2);
eventtWin(1,:) = [-1 2];   % for analysis of all trials
f = logspace(0,2.7,50);
numFreqs = length(f);
scalogram_metadata.f = f;
scalogram_metadata.precision = 'double';
filterBanks = cell(numTrialTypes, 1);

buffered_twin = eventtWin(1,:) + [-twin_buffer, twin_buffer];
numBufferedSamples = round(range(buffered_twin) * Fs);
filterBanks{1} = scalogramFilterBank(f, Fs, numBufferedSamples);
for iTrialType = 2 : numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
    buffered_twin = eventtWin(iTrialType,:) + [-twin_buffer, twin_buffer];
    numBufferedSamples = round(range(buffered_twin) * Fs);
    filterBanks{iTrialType} = scalogramFilterBank(f, Fs, numBufferedSamples);
end

for i_chDB = 1:length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    
    subject_scalogramDir = fullfile(scalogramDir, [implantID '_ps']);
    if ~exist(subject_scalogramDir,'dir')
        mkdir(subject_scalogramDir);
    end
    
%     lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
    lfp_directory = fullfile(lfp_root, [implantID '_HF_LFPs']);
    if ~exist(lfp_directory,'dir'); continue; end
    
    cd(lfp_directory);
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    if i_chDB==2
        startSession=27;
    else
        startSession = 1;
    end
    for iSession = startSession : numSessions
        session_scalogramDir = fullfile(subject_scalogramDir,[sessionList{iSession} '_scalograms']);
        if ~exist(session_scalogramDir,'dir')
            mkdir(session_scalogramDir);
        end
        fprintf('session %s, %d of %d\n', ...
            sessionList{iSession}, iSession, numSessions)
        scalogram_metadata.session = sessionList{iSession};
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.locationSubClass = ROI_list;
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        if isempty(sessionChannels);continue;end
        numSessionChannels = length(sessionChannels);
        
%         sessionComplete = true;
%         for iCh = 1 : numSessionChannels
%             for iTrialType = 1 : numTrialTypes
%                 trialType = trialTypeList{iTrialType};
%                 ch = sessionChannels{iCh};
%                 ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
%                 if ~exist(ch_scalogramName,'file');
%                     sessionComplete = false;
%                 end
%             end
%         end
%         if sessionComplete;continue;end
        
        % load the header for the current session
        sessionDate = sessionChannels{1}.date;
        if length(sessionDate) > 8
            sessionDate = datestr(sessionDate, 'yyyymmdd');
        end
        
        lfp_fileinfo = dir(['*_' sessionDate '*.hsdf']);
        if isempty(lfp_fileinfo); continue; end
        
        spikeInterp_idx = []; orig_lfp_idx = [];
        if length(lfp_fileinfo) > 1
            for iFile = 1 : length(lfp_fileinfo)
                if isempty(strfind(lfp_fileinfo(iFile).name, 'spikeInterp'))    % use only original files, ignore the "spike interpolation" files
                    orig_lfp_idx = [orig_lfp_idx, iFile];
                else
                    spikeInterp_idx = [spikeInterp_idx, iFile];
                end
            end
            lfp_fileinfo = lfp_fileinfo(orig_lfp_idx);
            if length(lfp_fileinfo) > 1
                disp([num2str(length(lfp_fileinfo)) ' files found for sessions on ' sessionDate]);
                break;
            end
        end
        
        lfp_fileName = lfp_fileinfo.name;
        
        % read in the LFP header and data
        header = getHSDHeader( lfp_fileName );
        numRecordingSites = header.main.num_channels;
        lfp_wireNums = zeros(1, numRecordingSites);
        for i_site = 1 : numRecordingSites
            lfp_wireNums(i_site) = header.channel(i_site).original_number;
        end
        
%         Fs = lfpFs( header );
%         startSample = round(Fs * twin_buffer);
        lfpDuration = getHSDlength( 'filename', lfp_fileName );
        
        disp(['loading ' lfp_fileName '...']);
%         tic
        lfp = readHSD( lfp_fileName, ...
                       numRecordingSites, ...
                       header.dataOffset, ...
                       Fs, ...
                       [0, lfpDuration] );
%         toc
        if isempty(lfp); continue; end
        
        trList = cell(numTrialTypes,1);
        for iTrialType = 1 : numTrialTypes
            trialType = trialTypeList{iTrialType};
            eventList = eventLists{iTrialType};
            twin = eventtWin(iTrialType,:);
            
            numSamples = round(range(twin) * Fs);
            
            ch = sessionChannels{1};    % just to get the trial list correct
            trialEventParams = getTrialEventParams(trialType);
            temp_trList = extractTrials(ch.trials, trialEventParams);

            % figure out how many VALID trials there are (that is, trials that don't
            % occur to close to the beginning or end of the session so that we can't
            % extract a full window around all events
            numTrials = length(temp_trList);
            minValidTrials = numTrials;
            numEvents = length(eventList);
            for iEventType = 1 : numEvents
                eventName = eventList{iEventType};
                
                numValidTrials = 1;
                for iTr = 1 : numTrials
                    if ~isfield(ch.trials(temp_trList(numValidTrials)).timestamps, eventName)
                        if numValidTrials == 1
                            temp_trList = temp_trList(2:end,:);
                        elseif numValidTrials == length(temp_trList)
                            temp_trList = temp_trList(1:numValidTrials);
                        else
                            temp_trList = [temp_trList(1:numValidTrials-1); temp_trList(numValidTrials+1:end)];
                        end
                        continue
                    end
                    event_ts = ch.trials(temp_trList(numValidTrials)).timestamps.(eventName);

                    tlim = event_ts + twin + [-twin_buffer, +twin_buffer];
                    % make sure the time window doesn't extend before the start of the
                    % recording or after the end of the recording
                    if any(tlim < 0) || any(tlim > lfpDuration)
                        if numValidTrials == 1
                            temp_trList = temp_trList(2:end,:);
                        elseif numValidTrials == length(temp_trList)
                            temp_trList = temp_trList(1:numValidTrials);
                        else
                            temp_trList = [temp_trList(1:numValidTrials-1); temp_trList(numValidTrials+1:end)];
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
            if numValidTrials == 0; continue; end
            trList{iTrialType} = temp_trList(1:numValidTrials);
        end
            
        for iCh = 1 : numSessionChannels
            ch = sessionChannels{iCh};

            if any(ch.wire.markedGood) == 0; continue; end

            repWire = getRepWire( sessionChannels{iCh} );

            lfp_idx = find( lfp_wireNums == repWire );
            if isempty(lfp_idx);continue;end

            volt_lfp = int2volt(lfp(lfp_idx, :), 'gain', header.channel(lfp_idx).gain);
            scalogram_metadata.channel = ch.name;
            
            fprintf('%s\n', ch.name);
            for iTrialType = 1 : numTrialTypes
                trialType = trialTypeList{iTrialType};
                numTrials = length(trList{iTrialType});
                if numTrials == 0; continue; end
                
                ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
%                 ch_scalogramName_md = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms_md.mat']);
                if exist(ch_scalogramName,'file');continue;end
                
                eventList = eventLists{iTrialType};
                twin = eventtWin(iTrialType, :);
                numSamples = round(range(twin) * Fs);
                buffered_twin = twin + [-twin_buffer, twin_buffer];
                
                scalogram_metadata.trialType = trialType;
                scalogram_metadata.eventList = eventList;
                scalogram_metadata.twin = twin;
                scalogram_metadata.t = linspace(twin(1),twin(2),numSamples);
                W = zeros(length(eventList),numSamples,numTrials,numFreqs);
                for iEvent = 1 : length(eventList)
                    lfp_events = extractLFP(volt_lfp, ch, trList{iTrialType}, eventList{iEvent}, buffered_twin, Fs);
%                     tic
%                     [full_W,~] = calculateComplexScalograms_EnMasse(lfp_events','Fs',Fs,'freqlist',f);
%                     toc
                    [full_W,~] = calculateComplexScalograms_EnMasse(lfp_events','filterbank',filterBanks{iTrialType});
                    W(iEvent,:,:,:) = full_W(startSample:startSample + numSamples - 1, :,:);
                    
                end    % for iEvent...
                
                save(ch_scalogramName, 'W', 'scalogram_metadata','-v7.3');
%                 save(ch_scalogramName_md, 'scalogram_metadata');
%                 fid = fopen(ch_scalogramName, 'w', machineFormat);
%                 fwrite(fid, real(W), scalogram_metadata.precision);
%                 fwrite(fid, imag(W), scalogram_metadata.precision);
%                 fclose(fid);
                
            end    % for iTrialType
        end    % for iCh...
        
    end    % for iSession...
    
end    % for i_chDB...
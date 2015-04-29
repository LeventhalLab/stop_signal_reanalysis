% script to calculate mean power spectra across trials for each
% tetrode-session

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
powerSpectraDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectra';
lfp_root          = ['/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs'];
sample_lfp_name = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs/IM164_LFPs/IM164_20091112_13-56-32.hsdf';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

twin = [-1 2];

eventList = {'noseCenterIn'};

ps_metadata.eventList = eventList;
ps_metadata.twin = twin;

for i_chDB = 5 : length(chDB_list)
    
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
    
    subject_powerSpectraDir = fullfile(powerSpectraDir, [implantID '_ps']);
    if ~exist(subject_powerSpectraDir,'dir')
        mkdir(subject_powerSpectraDir);
    end
    
    lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
    cd(lfp_directory);
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iSession = 1 : numSessions
        session_powerSpectraName = fullfile(subject_powerSpectraDir,['ps_' sessionList{iSession} '.mat']);
        if exist(session_powerSpectraName,'file')
            continue;
        end
        fprintf('session %s, %d of %d\n', ...
            sessionList{iSession}, iSession, numSessions)
        ps_metadata.session = sessionList{iSession};
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        cp = initChanParams();
        cp.tetrode = {'e2', 'e3'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numCh = length(sessionChannels);
        
        % load the header for the current session
        sessionDate = sessionChannels{1}.date;
        if length(sessionDate) > 8
            sessionDate = datestr(sessionDate, 'yyyymmdd');
        end
        
        lfp_fileinfo = dir(['*_' sessionDate '*.hsdf']);
        if isempty(lfp_fileinfo); continue; end
        
%         if length(lfp_fileinfo) ~= 1
%             disp([num2str(length(lfp_fileinfo)) ' files found for sessions on ' sessionDate]);
%             break;
%         end
        
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
        
        Fs = lfpFs( header );
        ps_metadata.Fs = Fs;
        
        lfpDuration = getHSDlength( 'filename', lfp_fileName );
        
            ch = sessionChannels{1};    % just to get the trial list correct
            trialEventParams = getTrialEventParams('any');
            trList = extractTrials(ch.trials, trialEventParams);

            % figure out how many VALID trials there are (that is, trials that don't
            % occur to close to the beginning or end of the session so that we can't
            % extract a full window around all events
            numTrials = length(trList);
            minValidTrials = numTrials;
            numEvents = length(eventList);
            for iEventType = 1 : numEvents
                eventName = eventList{iEventType};
                
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

                    tlim = event_ts + twin;
                    % make sure the time window doesn't extend before the start of the
                    % recording or after the end of the recording
                    if any(tlim < 0) || any(tlim > lfpDuration)
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
            if numValidTrials == 0; continue; end
            
        disp(['loading ' lfp_fileName '...']);
        tic
        lfp = readHSD( lfp_fileName, ...
                       numRecordingSites, ...
                       header.dataOffset, ...
                       Fs, ...
                       [0, lfpDuration] );
        toc
        if isempty(lfp); continue; end
        
        event_t = linspace(twin(1),twin(2),round(range(twin) * Fs));

        num_calculated_channels = 0;
        num_f = max(256,2^nextpow2(length(event_t)))/2+1;
        pxx = zeros(numCh,num_f,numValidTrials);
        ps_metadata.region = {};
        ps_metadata.chList = {};
        for iCh = 1:numCh%metadata.numWrittenChannels + 1 : numCh
%             iCh

            ch = sessionChannels{iCh};
            if any(ch.wire.markedGood) == 0; continue; end
            
            repWire = getRepWire( sessionChannels{iCh} );
            
            lfp_idx = find( lfp_wireNums == repWire );
            if isempty(lfp_idx);continue;end
            
            num_calculated_channels = num_calculated_channels + 1;
            ps_metadata.chList{iCh} = ch.name;
            ps_metadata.region{iCh} = ch.location.subclass;
            
            volt_lfp = int2volt(lfp(lfp_idx, :), 'gain', header.channel(lfp_idx).gain);
            
            lfp_events = extractLFP(volt_lfp, ch, trList, eventList{1}, twin, Fs);
            [pxx(iCh,:,:),f] = periodogram(lfp_events',[],[],Fs);
            
        end
        
        save(session_powerSpectraName,'pxx','f','ps_metadata');
        
    end
    
end
            
            
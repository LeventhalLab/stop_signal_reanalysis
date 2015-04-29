% script to check hilbert based power spectrograms
%%
chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
powerSpectrogramDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectrograms';
% phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';
lfp_root          = ['/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs'];
sample_lfp_name = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs/IM164_LFPs/IM164_20091112_13-56-32.hsdf';
[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventList{1} = {'noseCenterIn'};
eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideOut'};
eventList{3} = eventList{2};
eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

numTrialTypes = length(trialTypeList);

eventtWin(1,:) = [-1 2];   % for analysis of all trials
% analysisWin(1) = 3;
% stepSize(1)    = 3;
for iTrialType = 2 : numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
%     analysisWin(iTrialType) = 0.1;
%     stepSize(iTrialType)    = 0.05;
end

header = getHSDHeader( sample_lfp_name );
Fs = lfpFs(header);
M = [0.3 1 0.3];
test_freq = 20.5;
max_bw = 10;
freqBand = zeros(max_bw,4);
b = cell(max_bw,1);
test_freq = 19.5;
for i_bandwidth = 1 : max_bw
    freqBand(i_bandwidth,:) = test_freq + [-i_bandwidth-1.0 -i_bandwidth +i_bandwidth i_bandwidth+1.0];
    [n,fo,mo,w] = firpmord(freqBand(i_bandwidth,:), ...
                                           M, ...
                                           [1 1 1] * 0.1, ...
                                           Fs);
    b{i_bandwidth} = firpm(n,fo,mo,w);
end


%%
iTrialType = 2;
twin = eventtWin(iTrialType,:);
trialType = trialTypeList{iTrialType};
rowSpace = 0.5;
for i_chDB = 2 : 2%length(chDB_list)
    
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
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    
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
    
    for iSession = numSessions-1 : numSessions
        
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
        lfpDuration = getHSDlength( 'filename', lfp_fileName );
        
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
            
            % find the first index of the 1 Hz bins
            startIdx_1Hz = find((centerFreqs_1Hz > max025Hz), 1,'first');
            
            f = [centerFreqs_025Hz; centerFreqs_1Hz(startIdx_1Hz:end)];

            num_freq_025 = length(centerFreqs_025Hz);

            Fs = md_025Hz.metadata.Fs;
            powerSpect_metadata.Fs = Fs;
            numFreqs = length(f);
        
            
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
        
        disp(['loading ' lfp_fileName '...']);
        tic
        lfp = readHSD( lfp_fileName, ...
                       numRecordingSites, ...
                       header.dataOffset, ...
                       Fs, ...
                       [0, lfpDuration] );
        toc
        if isempty(lfp); continue; end
        
        % read in the LFP header and data
        header = getHSDHeader( lfp_fileName );
        numRecordingSites = header.main.num_channels;
        lfp_wireNums = zeros(1, numRecordingSites);
        for i_site = 1 : numRecordingSites
            lfp_wireNums(i_site) = header.channel(i_site).original_number;
        end
        
        Fs = lfpFs( header );
        lfpDuration = getHSDlength( 'filename', lfp_fileName );
        
        % create a text file with metadata for the files that will be written
        % to this folder
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir_1Hz, metadata_filename);
        
        if exist(metadata_filename, 'file')
            md_1Hz = load(metadata_filename);
%         else
%             metadata.session   = sessionList{iSession};
%             metadata.duration  = lfpDuration;
%             metadata.numSamps  = 0;
%             metadata.sigmax    = zeros(1, numCh); 
%             metadata.Fs        = Fs;
%             metadata.chNames   = cell(1, numCh);
%             metadata.freqBands = freqBands;
%             metadata.bitOrder  = bitOrder;
%             for iCh = 1 : numCh
%                 metadata.chNames{iCh} = sessionChannels{iCh}.name;
%             end
% 
%             metadata.numWrittenChannels = 0;   % number of channels for which all frequencies have been filtered, hilbert transformed, and stored
%             metadata.numWrittenFreqs    = 0;
% 
%             metadata.b = cell(1, numBands);
%             metadata.comment = ['Analytic signals calculated by bandpass filtering LFPs and calculating the Hilbert transform. ' ...
%                                 'Data are stored in binary files as double precision complex numbers. Each analytic signal is ' ...
%                                 'numSamps values long. Analytic signals are stored column-wise, with the left-most column representing ' ...
%                                 'the LFP filtered in the freqBands(:,1) band, etc. The complex numbers are stored as 2 sets of double ' ...
%                                 'precision numbers - first, the full stream of real parts, then the full stream of imaginary parts '...
%                                 '(that is, [1+2i, 3+4i] would be stored as [1,3,2,4]. Double precision numbers were scaled to the full ' ...
%                                 'int16 range and written as int16s. That is if x is a 16-bit integer read from the file, to get the '...
%                                 'correct value, take x * sigmax(channel_index) / 32767. sigmax is an array stored in the metadata ' ...
%                                 'structure.'];
%             save(metadata_filename, 'metadata');
        end
        
        numSamples = size(lfp, 2);

        t = linspace(1/Fs, lfpDuration, numSamples);
        event_t = linspace(twin(1),twin(2),round(range(twin) * Fs));
        chtic = tic;
        
            ch = sessionChannels{1};    % just to get the trial list correct
            trialEventParams = getTrialEventParams(trialType);
            trList = extractTrials(ch.trials, trialEventParams);

            % figure out how many VALID trials there are (that is, trials that don't
            % occur to close to the beginning or end of the session so that we can't
            % extract a full window around all events
            numTrials = length(trList);
            minValidTrials = numTrials;
            numEvents = length(eventList{iTrialType});
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
            
            
        for iCh = 3:numCh%metadata.numWrittenChannels + 1 : numCh
            iCh
            disp(sprintf('trialtype: %s, session %s, %d of %d; channel %d of %d', ...
                trialType, sessionList{iSession}, iSession, numSessions, iCh, numCh))

            ch = sessionChannels{iCh};
            if any(ch.wire.markedGood) == 0; continue; end
            
            repWire = getRepWire( sessionChannels{iCh} );
            
            lfp_idx = find( lfp_wireNums == repWire );
            if isempty(lfp_idx);continue;end
            
            volt_lfp = int2volt(lfp(lfp_idx, :), 'gain', header.channel(lfp_idx).gain);
            filt_lfp = zeros(max_bw,length(volt_lfp));
            for i_bandwidth = 1 : max_bw
                filt_lfp(i_bandwidth,:) = filtfilt(b{i_bandwidth},1,volt_lfp);
            end
            for iEvent = 4 : numEvents
%                 tic
                for iFreq = 46 : numFreqs
% iFreq
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
                    powerSpect(iEvent, iFreq, :) = squeeze(mean(freqPower,1));
                    
                    lfp_events = extractLFP(volt_lfp, ch, trList, eventList{iTrialType}{iEvent}, twin, Fs);
                    filt_events = zeros(max_bw,size(lfp_events,1),size(lfp_events,2));
                    for i_bandwidth = 1 : max_bw
                        filt_events(i_bandwidth,:,:) = extractLFP(filt_lfp(i_bandwidth,:), ch, trList, eventList{iTrialType}{iEvent}, twin, Fs);
                    end
                    
                    figure
%                     for ii = 1 : size(lfp_events,1)
%                         plot(event_t,lfp_events(ii,:) + ii*rowSpace,'color','b');
%                         hold on
%                         plot(event_t,abs(ansig(ii,:)) + ii*rowSpace,'color','r');
%                         plot(event_t,real(ansig(ii,:)) + ii*rowSpace,'color','g');
%                     end
                    plot(event_t,lfp_events(end,:),'color','b');
                    hold on
%                     for i_bandwidth = 1 : 4 : 5%max_bw
%                         plot(event_t,squeeze(filt_events(i_bandwidth,end,:)));
%                     end
                    plot(event_t,squeeze(filt_events(1,end,:)),'color','r');
                    plot(event_t,squeeze(filt_events(5,end,:)),'color','k');
                    legend('raw signal','1 Hz bandwidth','5 Hz bandwidth');
                    title('center frequency 20.5 Hz')
                end
%                 toc
            end
            
        end
        
    end
    
end
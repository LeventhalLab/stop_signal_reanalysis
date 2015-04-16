% script_calc_power_RTcorrelations

% UPDATED 10-27-2014

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
% phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};
numEvents = length(eventList);
twin = [-1 1];

trialType = 'correctgo';

for i_chDB = 4 : 4%length(chDB_list)

    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir, 'dir')
        disp([subject_hilbertDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    subject_powerRTcorrdir = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr']);
    if ~exist(subject_powerRTcorrdir, 'dir')
        mkdir(subject_powerRTcorrdir);
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
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
    allRT = sort(allRT);
    numRT = length(allRT);
    
    powerRTcorr = cell(1, numSessions);    % cell array to store correlations between continuous narrow-band power and RT

    for iSession = 26:26%1 : numSessions

        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        numCh = length(sessionChannels);

        hilbert_sessionDir   = fullfile(subject_hilbertDir, sessionList{iSession});
        if ~exist(hilbert_sessionDir, 'dir'); continue; end
        
        powerRTcorr_sessionDir = fullfile(subject_powerRTcorrdir, sessionList{iSession});
        if ~exist(powerRTcorr_sessionDir, 'dir')
            mkdir(powerRTcorr_sessionDir);
        end
        
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
        
        if ~exist(metadata_filename, 'file')
            error([metadata_filename ' could not be found.']);
        end
        load(metadata_filename);
        
        centerFreqs = mean(metadata.freqBands, 2);
        numFreqs = length(centerFreqs);
        
        powerRTcorr_metadata.freqs = centerFreqs;
        powerRTcorr_metadata.chNames = metadata.chNames;
        powerRTcorr_metadata.eventList = eventList;
        powerRTcorr_metadata.twin = twin;
        powerRTcorr_metadata.Fs = metadata.Fs;
        powerRTcorr_metadata.trialType = trialType;
        powerRTcorr_metadata.completeChannels = 0;

        numSamps = round(range(twin) * metadata.Fs);
        
        for iCh = 1 : numCh
            disp(sprintf('trialtype: %s, session %s, %d of %d; channel %d of %d', ...
                trialType, sessionList{iSession}, iSession, numSessions, iCh, numCh))
            ch = sessionChannels{iCh};
            
            hilbert_name = ['analytic_' sessionChannels{iCh}.name '.bin'];
            hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
            
            power_RTcorr_name = ['power_RTcorr_' sessionChannels{iCh}.name '.mat'];
            power_RTcorr_name = fullfile(powerRTcorr_sessionDir, power_RTcorr_name);
            
            powerRTcorr = NaN(numEvents, numFreqs, numSamps);
            
            % load analytic signal around each event and calculate
            % correlation coefficients
            if exist(power_RTcorr_name, 'file')
                continue
            end
            
            for iEvent = 1 : numEvents
                tic
                for iFreq = 1 : numFreqs
                    
                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    iFreq, ...
                                                    eventList{iEvent}, ...
                                                    trialType, ...
                                                    twin, ...
                                                    'hilbertdir', hilbert_directory);
                    freqPower = abs(ansig) .^ 2;

                    for iSamp = 1 : numSamps
                        [cc, p] = corr(freqPower(:, iSamp), RT{iSession}', ...
                                       'type','spearman');
                        powerRTcorr(iEvent, iFreq, iSamp) = cc;
                    end
                        
                end
                toc
            end

            save(power_RTcorr_name, 'powerRTcorr', 'powerRTcorr_metadata');
            
        end
    end
    
end
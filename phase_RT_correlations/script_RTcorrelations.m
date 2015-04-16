% script_RTcorrelations

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};
numEvents = length(eventList);
twin = [-1 1];

% parameters for phase-RT analysis; adapted from VanRullen et al, Ongoing
% EEG phase..., Frontiers in Psychology, 2011; details in Drewes and
% Vanrullen, "This is the rhtyhm of your eyes: the phase of ongoing
% electroencephalogram oscillations modulates saccadic reaction time.", J
% Neurosci, 2011
RTquantiles = 0.2 : 0.2 : 1;    % quintiles for dividing up the RT distribution
numQuantiles = length(RTquantiles);
trialType = 'correctgo';

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir, 'dir')
        disp([subject_hilbertDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    subject_powerRTcorrdir = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr']);
    if ~exist(subject_powerRTcorrdir, 'dir')
        mkdir(subject_powerRTcorrdir);
    end
    subject_phaseRTcorrdir = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr']);
    if ~exist(subject_phaseRTcorrdir, 'dir')
        mkdir(subject_phaseRTcorrdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
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
    quant_idx = round(RTquantiles * numRT);
    RTquantile_borders = [0, allRT(quant_idx)];
    
    powerRTcorr = cell(1, numSessions);    % cell array to store correlations between continuous narrow-band power and RT
    RTquantile_phases = cell(1, numSessions);    % cell array to store phases across RT quantiles
    
    for iSession = 1 : numSessions

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
        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorrdir, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            mkdir(phaseRTcorr_sessionDir);
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
        
        phaseRTcorr_metadata.freqs = centerFreqs;
        phaseRTcorr_metadata.chNames = metadata.chNames;
        phaseRTcorr_metadata.eventList = eventList;
        phaseRTcorr_metadata.twin = twin;
        phaseRTcorr_metadata.Fs = metadata.Fs;
        phaseRTcorr_metadata.trialType = trialType;
        phaseRTcorr_metadata.RTquantiles = RTquantiles;
        phaseRTcorr_metadata.RTcutoffs = RTquantile_borders;
        
        numSamps = round(range(twin) * metadata.Fs);
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            hilbert_name = ['analytic_' sessionChannels{iCh}.name '.bin'];
            hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
            
            phase_RTcorr_name = ['phase_RT_analysis_' sessionChannels{iCh}.name '.mat'];
            phase_RTcorr_name = fullfile(phaseRTcorr_sessionDir, phase_RTcorr_name);
            
            power_RTcorr_name = ['power_RTcorr_' sessionChannels{iCh}.name '.mat'];
            power_RTcorr_name = fullfile(powerRTcorr_sessionDir, power_RTcorr_name);
            
            powerRTcorr = NaN(numEvents, numFreqs, numSamps);
            RTphases    = cell(numEvents, numFreqs, numQuantiles);
            
            % load analytic signal around each event and calculate
            % correlation coefficients
            if exist(phase_RTcorr_name, 'file') && exist(power_RTcorr_name, 'file')
                continue
            end
            
            for iEvent = 1 : numEvents
                tic
                for iFreq = 1 : numFreqs
                    
                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    iFreq, ...
                                                    eventList{iEvent}, ...
                                                    trialType, ...
                                                    twin );
                    freqPower = abs(ansig) .^ 2;
                    freqPhase = angle(ansig);

                    for iQuant = 1 : numQuantiles
                        quantIdx = (RT{iSession} > RTquantile_borders(iQuant) & ...
                                    RT{iSession} < RTquantile_borders(iQuant + 1));
                        RTphases{iEvent, iFreq, iQuant} = freqPhase(quantIdx, :);
                    end

                    for iSamp = 1 : numSamps
                        [cc, p] = corr(freqPower(:, iSamp), RT{iSession}', ...
                                       'type','spearman');
                        powerRTcorr(iEvent, iFreq, iSamp) = cc;
                    end
                        
                end
                toc
            end
            
            save(phase_RTcorr_name, 'RTphases', 'phaseRTcorr_metadata');
            save(power_RTcorr_name, 'powerRTcorr', 'powerRTcorr_metadata');
            
        end
    end
    
end
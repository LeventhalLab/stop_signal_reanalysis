% script_analyzeRTphasehists

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

num_phaseBins = 12;
phaseBins = linspace(-pi+pi/num_phaseBins, pi-pi/num_phaseBins, num_phaseBins);
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

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
    
    subject_phaseRThistdir = fullfile(phaseRThist_directory, [implantID '_phaseRTcorr']);
    if ~exist(subject_phaseRThistdir, 'dir')
        continue;
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
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
        
        phaseRThist_sessionDir = fullfile(subject_phaseRThistdir, sessionList{iSession});
        if ~exist(phaseRThist_sessionDir, 'dir')
            continue;
        end
        
        for iCh = 1 : numCh
            iCh
            ch = sessionChannels{iCh};
            saveName = ['phase_RT_circCorr_' ch.name '.mat'];
            saveName = fullfile(phaseRThist_sessionDir, saveName);
            if exist(saveName, 'file'); continue; end
            
            phase_RThist_name = ['phase_RT_hist_analysis_' ch.name '.mat'];
            phase_RThist_name = fullfile(phaseRThist_sessionDir, phase_RThist_name);
            if ~exist(phase_RThist_name, 'file'); continue; end
            
            load(phase_RThist_name);
            numEvents = length(phaseRThist_metadata.eventList);
            numFreqs = length(phaseRThist_metadata.freqs);
            numSamps = size(freqPhase, 4);
            pcc = zeros(numEvents, numSamps, numFreqs);
            for iEvent = 1 : numEvents
                iEvent
                for iFreq = 1 : numFreqs
                    for i_t = 1 : numSamps  
                        % calculate the circular correlation coefficient
                        % between phase and RT
                        phaseVec = squeeze(freqPhase(iEvent, iFreq, :, i_t));
                        pcc(iEvent, i_t, iFreq) = circ_corrcl(phaseVec, RT);  % circular-linear pearson correlation coefficient
                    end
                end
            end
            save(saveName, 'pcc', 'phaseRThist_metadata');
        end    % for iCh = 1 : numCh

    end
    
end
                        
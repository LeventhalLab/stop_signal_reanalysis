chDB_directory    = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
phase_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_power_gabors';

numSurrogates = 200;

for i_chDB = 1 : 4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
        
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        disp([subject_stopPhaseDir ' not found. Skipping ' implantID '...'])
        continue
    end

    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    cp = initChanParams();
    cp.task = 3;
    chList = extractChannels(cp, channels);
    channels = channels(chList);
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iSession = 1 : numSessions
        
        disp(sessionList{iSession})
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.task = 3;    % only stop-signal sessions
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        cp = initChanParams();
        cp.tetrode = {'e2', 'e3', 'e02','e03'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numSessionChannels = length(sessionChannels);
        
        regionList = getRegionsfromChannelDB(sessionChannels);
        numRegions = length(regionList);
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end
        
        for iCh = 1 : numSessionChannels
            
                ch = sessionChannels{iCh};
                surr_mrv_mat_saveName = ['mrv_stopSuccvFail_surrogates_' ch.name '_gabor.mat'];
                surr_mrv_mat_saveName = fullfile(session_stopPhaseDir, surr_mrv_mat_saveName);
                if ~exist(surr_mrv_mat_saveName, 'file')
                    continue
                end
                
                load(surr_mrv_mat_saveName);
                surr_STOPmrv_metadata.numSurrogates = numSurrogates;
                
                save(surr_mrv_mat_saveName, 'mean_surr_diff','std_surr_diff','surr_STOPmrv_metadata');
                
        end
        
    end
    
end
                
                
                
                
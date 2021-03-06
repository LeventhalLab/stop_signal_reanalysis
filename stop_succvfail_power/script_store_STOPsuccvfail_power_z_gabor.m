% script_store_STOPsuccvfail_power_z_gabor

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

numSurrogates = 200;

for i_chDB = 1:4%4%length(chDB_list)
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_stopPowerDir = fullfile(power_stop_succvfail_directory, [implantID '_stopPower']);
    if ~exist(subject_stopPowerDir, 'dir')
        disp([subject_stopPowerDir ' not found. Skipping ' implantID '...'])
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
        
%         subject_stopPowerDir = fullfile(power_stop_succvfail_directory, [implantID '_stopPower']);
        
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
        
        session_stopPowerDir = fullfile(subject_stopPowerDir, sessionList{iSession});
        if ~exist(session_stopPowerDir, 'dir')
            disp([session_stopPowerDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end
        
        for iCh = 1 : numSessionChannels
            stopPower_name = ['stopPower_' sessionChannels{iCh}.name '_gabor.mat'];
            stopPower_name = fullfile(session_stopPowerDir, stopPower_name);
            if ~exist(stopPower_name, 'file'); continue; end
            
            load(stopPower_name);
            break;
        end
        
        if ~exist(stopPower_name, 'file'); continue; end
        z_STOPpower_metadata.Fs = STOPmetadata.Fs;
        z_STOPpower_metadata.f = STOPmetadata.f;
        z_STOPpower_metadata.eventList = STOPmetadata.eventList;
        z_STOPpower_metadata.chName = STOPmetadata.chName;
        z_STOPpower_metadata.twin = STOPmetadata.twin;
        
        numEventTypes = length(STOPmetadata.eventList);
        numFreqs  = length(STOPmetadata.f);
        numSamps  = size(correctSTOP_power, 2);
        
        powerDiff = zeros(numEventTypes, numSamps, numFreqs);
        t = linspace(z_STOPpower_metadata.twin(1),z_STOPpower_metadata.twin(2),numSamps);
        f = STOPmetadata.f;
        
        z_STOPpower_metadata.t = t;
        z_STOPpower_metadata.f = f;
        
        mean_power = zeros(numRegions, 2, numEventTypes, numSamps, numFreqs);   % 1 for correct, 2 for failed STOP
        mean_powerDiff = zeros(numRegions, numEventTypes, numSamps, numFreqs);  
        num_ch_per_region = zeros(1, numRegions);
        
        for iRegion = 1 : numRegions
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numRegionChannels = length(regionChannels);
            
            power_by_region{iRegion} = zeros(numRegionChannels, 2, numEventTypes, numSamps, numFreqs);
            
            num_ch_per_region(iRegion) = numRegionChannels;

            for iCh = 1 : numRegionChannels
%                 iCh
            
                ch = regionChannels{iCh};
                z_STOPpower_metadata.chName = ch.name;
                z_STOPpower_metadata.region = ch.location.name;
                surr_power_mat_saveName = ['power_stopSuccvFail_surrogates_' ch.name '_gabor.mat'];
                surr_power_mat_saveName = fullfile(session_stopPowerDir, surr_power_mat_saveName);
                z_power_mat_saveName = ['power_stopSuccvFail_z_' ch.name '_gabor.mat'];
                z_power_mat_saveName = fullfile(session_stopPowerDir, z_power_mat_saveName);
                if exist(z_power_mat_saveName, 'file')
                    continue
                end
                if ~exist(surr_power_mat_saveName, 'file')
                    continue
                end
                load(surr_power_mat_saveName);
                
                stopPower_name = ['stopPower_' ch.name '_gabor.mat'];
                stopPower_name = fullfile(session_stopPowerDir, stopPower_name);
                if ~exist(stopPower_name, 'file'); continue; end
                load(stopPower_name);
                
                numEventTypes = size(correctSTOP_power, 1);
                numFreqs      = size(correctSTOP_power, 4);
                numSamps      = size(correctSTOP_power, 2);
                numSuccTrials = size(correctSTOP_power, 3);
                numFailTrials = size(failedSTOP_power, 3);
                totSTOPTrials = numSuccTrials + numFailTrials;
                
                realPowerDiff = squeeze(mean(correctSTOP_power, 3)) - squeeze(mean(failedSTOP_power,3));
                zPowerDiff = (realPowerDiff - mean_surr_diff) ./ std_surr_diff;
                
                save(z_power_mat_saveName, 'realPowerDiff','zPowerDiff','z_STOPpower_metadata');
                
            end
            
        end
        
    end
    
end
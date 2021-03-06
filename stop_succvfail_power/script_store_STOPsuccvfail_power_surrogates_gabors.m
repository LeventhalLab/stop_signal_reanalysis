% script_store_STOPsuccvfail_power_surrogates_gabors

% script to go through all STOP sessions and generate surrogate
% distributions for differences in the mrv's between STOP-success and
% STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

numSurrogates = 200;

for i_chDB = 2 : 4%length(chDB_list)
    
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
        surr_STOPpower_metadata.Fs = STOPmetadata.Fs;
        surr_STOPpower_metadata.f = STOPmetadata.f;
        surr_STOPpower_metadata.eventList = STOPmetadata.eventList;
        surr_STOPpower_metadata.chName = STOPmetadata.chName;
        surr_STOPpower_metadata.twin = STOPmetadata.twin;
        
        numEventTypes = length(STOPmetadata.eventList);
        numFreqs  = length(STOPmetadata.f);
        numSamps  = size(correctSTOP_power, 2);

        powerDiff = zeros(numEventTypes, numFreqs, numSamps);
        t = linspace(surr_STOPpower_metadata.twin(1),surr_STOPpower_metadata.twin(2),numSamps);
        f = STOPmetadata.f;
        
        surr_STOPpower_metadata.t = t;
        
        surrogate_vecDiffmat_saveName = ['vecDiff_stopSuccvFail_surrogates_' sessionList{iSession} '_gabor.mat'];
        surrogate_vecDiffmat_saveName = fullfile(session_stopPowerDir, surrogate_vecDiffmat_saveName);

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
                surr_STOPpower_metadata.chName = ch.name;
                surr_STOPpower_metadata.region = ch.location.name;
                surr_power_mat_saveName = ['power_stopSuccvFail_surrogates_' ch.name '_gabor.mat'];
                surr_power_mat_saveName = fullfile(session_stopPowerDir, surr_power_mat_saveName);
                if exist(surr_power_mat_saveName, 'file')
                    continue
                end
            
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
                
                allPower = zeros(numEventTypes, numSamps, totSTOPTrials, numFreqs);
                for iEventType = 1 : numEventTypes
                    for iFreq = 1 : numFreqs
                        allPower(iEventType, :, :, iFreq) = [squeeze(correctSTOP_power(iEventType, :, :, iFreq)), ...
                                                             squeeze(failedSTOP_power(iEventType, :, :, iFreq))];
                    end
                end
                surr_diff = zeros(numSurrogates, numEventTypes, numSamps, numFreqs);
                for iSurrogate = 1 : numSurrogates
                    surrogateIndices = randperm(totSTOPTrials);
                    
                    surr_correctPower = allPower(:, :, surrogateIndices(1:numSuccTrials), :);
                    surr_failedPower = allPower(:, :, surrogateIndices(numSuccTrials+1:end), :);
                    
                    surr_power1 = squeeze(mean(surr_correctPower, 3));
                    surr_power2 = squeeze(mean(surr_failedPower, 3));

                    surr_diff(iSurrogate, :, :, :) = surr_power1 - surr_power2;
                end
                
                mean_surr_diff = squeeze(mean(surr_diff, 1));
                std_surr_diff  = squeeze(std(surr_diff, 0, 1));
                
                save(surr_power_mat_saveName, 'mean_surr_diff','std_surr_diff','surr_STOPpower_metadata');
                    
            end    % for iCh...
            
        end    % for iRegion...
        
    end    % for iSession...
    
end
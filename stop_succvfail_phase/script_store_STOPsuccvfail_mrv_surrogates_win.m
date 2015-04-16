% script_store_STOPsuccvfail_mrv_surrogates_b

% script to go through all STOP sessions and generate surrogate
% distributions for differences in the mrv's between STOP-success and
% STOP-failure trials

chDB_directory    = '\\172.20.138.142\PublicLeventhal1\dan\stop-signal reanalysis\stop-signal data structures';
% hilbert_1Hz_directory = '\\172.20.138.142\PublicLeventhal1\dan\stop-signal reanalysis\Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '\\172.20.138.142\PublicLeventhal1\dan\stop-signal reanalysis\Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '\\172.20.138.142\PublicLeventhal1\dan\stop-signal reanalysis\stop_succvfail_phase';
power_stop_succvfail_directory = '\\172.20.138.142\PublicLeventhal1\dan\stop-signal reanalysis\stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

numSurrogates = 200;

for i_chDB = 2 : 2%length(chDB_list)
    
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
        
        disp(sprintf('%s, %d of %d', sessionList{iSession}, iSession, numSessions))
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.task = 3;    % only stop-signal sessions
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
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
            stopPhase_name = ['stopPhases_' sessionChannels{iCh}.name '.mat'];
            stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
            if ~exist(stopPhase_name, 'file'); continue; end
            
            load(stopPhase_name);
            break;
        end
        if ~exist(stopPhase_name, 'file'); continue; end
        surr_STOPmrv_metadata.Fs = STOPmetadata.Fs;
        surr_STOPmrv_metadata.freqBands = STOPmetadata.freqBands;
        surr_STOPmrv_metadata.eventList = STOPmetadata.eventList;
        surr_STOPmrv_metadata.chName = STOPmetadata.chName;
        surr_STOPmrv_metadata.twin = STOPmetadata.twin;
        
        numEventTypes = length(STOPmetadata.eventList);
        numFreqs  = size(STOPmetadata.freqBands, 1);
        numSamps  = size(correctSTOP_phase, 4);

        mrvDiff = zeros(numEventTypes, numFreqs, numSamps);
        t = linspace(surr_STOPmrv_metadata.twin(1),surr_STOPmrv_metadata.twin(2),numSamps);
        f = mean(STOPmetadata.freqBands, 2);
        
        surr_STOPmrv_metadata.t = t;
        surr_STOPmrv_metadata.f = f;
%         surr_STOPmrv_metadata.regionList = regionList;
%         surr_STOPmrv_metadata.chNames = cell(1, numRegions);
        
        surrogate_vecDiffmat_saveName = ['vecDiff_stopSuccvFail_surrogates_' sessionList{iSession} '.mat'];
        surrogate_vecDiffmat_saveName = fullfile(session_stopPhaseDir, surrogate_vecDiffmat_saveName);
%         if exist(vecDiffmat_saveName,'file')
%             load(vecDiffmat_saveName);
%             if STOPmrv_metadata.num_regions_complete == numRegions; continue; end
%             num_ch_per_region = STOPmrv_metadata.num_ch_per_region;
%         else
%             mrvDiff_by_region = cell(1, numRegions);
            mean_mrvDiff = zeros(numRegions, numEventTypes, numFreqs, numSamps);
            mean_mrv = zeros(numRegions, 2, numEventTypes, numFreqs, numSamps);   % 1 for correct, 2 for failed STOP
            mean_mrvDiff = zeros(numRegions, numEventTypes, numFreqs, numSamps);  
            num_ch_per_region = zeros(1, numRegions);
%             STOPmrv_metadata.num_regions_complete = 0;
%         end

        for iRegion = 1 : numRegions
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numRegionChannels = length(regionChannels);
            
%             mrvDiff_by_region{iRegion} = zeros(numRegionChannels, numEventTypes, numFreqs, numSamps);
            mrv_by_region{iRegion} = zeros(numRegionChannels, 2, numEventTypes, numFreqs, numSamps);
            
            num_ch_per_region(iRegion) = numRegionChannels;
%             STOPmrv_metadata.chNames{iRegion} = cell(1, numRegionChannels);

            for iCh = 1 : numRegionChannels
            
                ch = regionChannels{iCh};
                surr_STOPmrv_metadata.chName = ch.name;
                surr_STOPmrv_metadata.region = ch.location.name;
                surr_mrv_mat_saveName = ['mrv_stopSuccvFail_surrogates_' ch.name '.mat'];
                surr_mrv_mat_saveName = fullfile(session_stopPhaseDir, surr_mrv_mat_saveName);
                if exist(surr_mrv_mat_saveName, 'file')
                    continue
                end
            
                stopPhase_name = ['stopPhases_' ch.name '.mat'];
                stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
                if ~exist(stopPhase_name, 'file'); continue; end

                load(stopPhase_name);
                numEventTypes = size(correctSTOP_phase, 1);
                numFreqs      = size(correctSTOP_phase, 2);
                numSamps      = size(correctSTOP_phase, 4);
                numSuccTrials = size(correctSTOP_phase, 3);
                numFailTrials = size(failedSTOP_phase, 3);
                totSTOPTrials = numSuccTrials + numFailTrials;
                
                allPhases = zeros(numEventTypes, numFreqs, totSTOPTrials, numSamps);
                for iEventType = 1 : numEventTypes
                    for iFreq = 1 : numFreqs
                        allPhases(iEventType, iFreq, :, :) = [squeeze(correctSTOP_phase(iEventType, iFreq, :, :)); ...
                                                              squeeze(failedSTOP_phase(iEventType, iFreq, :, :))];
                    end
                end
                surr_diff = zeros(numSurrogates, numEventTypes, numFreqs, numSamps);
                for iSurrogate = 1 : numSurrogates
                    surrogateIndices = randperm(totSTOPTrials);
                    
                    surr_correctPhase = allPhases(:, :, surrogateIndices(1:numSuccTrials), :);
                    surr_failedPhase = allPhases(:, :, surrogateIndices(numSuccTrials+1:end), :);

                    surr_mrv1 = squeeze(mean(exp(1i*surr_correctPhase), 3));
                    surr_mrv2 = squeeze(mean(exp(1i*surr_failedPhase), 3));
                    surr_diff(iSurrogate, :, :, :) = abs(surr_mrv1 - surr_mrv2);
                end
                
                mean_surr_diff = squeeze(mean(surr_diff, 1));
                std_surr_diff  = squeeze(std(surr_diff, 0, 1));
                
                save(surr_mrv_mat_saveName, 'mean_surr_diff','std_surr_diff','surr_STOPmrv_metadata');
                    
            end    % for iCh...
            
        end    % for iRegion...
        
    end    % for iSession...
    
end
% script_analyze_STOPsuccvfail_phases

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

for i_chDB = 1 : 1%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        disp([subject_stopPhaseDir ' not found. Skipping ' implantID '...'])
        continue
    end

    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iSession = 3 : 3%numSessions
        
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
        STOPmrv_metadata.Fs = STOPmetadata.Fs;
        STOPmrv_metadata.freqBands = STOPmetadata.freqBands;
        STOPmrv_metadata.eventList = STOPmetadata.eventList;
        STOPmrv_metadata.chName = STOPmetadata.chName;
        STOPmrv_metadata.twin = STOPmetadata.twin;
        
        numEventTypes = length(STOPmetadata.eventList);
        numFreqs  = size(STOPmetadata.freqBands, 1);
        numSamps  = size(correctSTOP_phase, 4);

        mrvDiff = zeros(numEventTypes, numFreqs, numSamps);
        t = linspace(STOPmrv_metadata.twin(1),STOPmrv_metadata.twin(2),numSamps);
        f = mean(STOPmetadata.freqBands, 2);
        
        STOPmrv_metadata.t = t;
        STOPmrv_metadata.f = f;
        STOPmrv_metadata.regionList = regionList;
        STOPmrv_metadata.chNames = cell(1, numRegions);
        
        vecDiffmat_saveName = ['vecDiff_stopSuccvFail_' sessionList{iSession} '.mat'];
        vecDiffmat_saveName = fullfile(session_stopPhaseDir, vecDiffmat_saveName);
%         if exist(vecDiffmat_saveName,'file')
%             load(vecDiffmat_saveName);
%             if STOPmrv_metadata.num_regions_complete == numRegions; continue; end
%             num_ch_per_region = STOPmrv_metadata.num_ch_per_region;
%         else
%             mrvDiff_by_region = cell(1, numRegions);
            mrv_by_region = cell(1, numRegions);
            mean_mrvDiff = zeros(numRegions, numEventTypes, numFreqs, numSamps);
            mean_mrv = zeros(numRegions, 2, numEventTypes, numFreqs, numSamps);   % 1 for correct, 2 for failed STOP
            mean_mrvDiff = zeros(numRegions, numEventTypes, numFreqs, numSamps);  
            num_ch_per_region = zeros(1, numRegions);
            STOPmrv_metadata.num_regions_complete = 0;
%         end

        for iRegion = STOPmrv_metadata.num_regions_complete + 1 : numRegions
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numRegionChannels = length(regionChannels);
            
%             mrvDiff_by_region{iRegion} = zeros(numRegionChannels, numEventTypes, numFreqs, numSamps);
            mrv_by_region{iRegion} = zeros(numRegionChannels, 2, numEventTypes, numFreqs, numSamps);
            
            num_ch_per_region(iRegion) = numRegionChannels;
            STOPmrv_metadata.chNames{iRegion} = cell(1, numRegionChannels);

            for iCh = 1 : numRegionChannels
                iCh
            
                ch = regionChannels{iCh};
                STOPmrv_metadata.chNames{iRegion}{iCh} = ch.name;
                mrv_mat_saveName = ['mrv_stopSuccvFail_' ch.name '.mat'];
                mrv_mat_saveName = fullfile(session_stopPhaseDir, mrv_mat_saveName);
                if exist(mrv_mat_saveName, 'file')
                    load(mrv_mat_saveName);
                    mrv_by_region{iRegion}(iCh, :, :, :, :) = mrv;
                    continue
                end
            
                stopPhase_name = ['stopPhases_' ch.name '.mat'];
                stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
                if ~exist(stopPhase_name, 'file'); continue; end

                load(stopPhase_name);
            
                mrv(1, :, :, :) = squeeze(mean(exp(1i*correctSTOP_phase), 3));
                mrv(2, :, :, :) = squeeze(mean(exp(1i*failedSTOP_phase), 3));
                mrv_by_region{iRegion}(iCh, :, :, :, :) = mrv;
                
                re_mrv = real(mrv); im_mrv = imag(mrv);
                
                save(mrv_mat_saveName, 're_mrv','im_mrv','STOPmrv_metadata');
                    
            end    % for iCh...
            mean_mrv(iRegion, :, :, :, :) = squeeze(mean(mrv_by_region{iRegion}, 1));
            re_mean_mrv = real(mean_mrv);
            im_mean_mrv = imag(mean_mrv);
%             mean_mrvDiff(iRegion, :, : ,:) = squeeze(mean(mrvDiff_by_region{iRegion}, 1));
            
            STOPmrv_metadata.num_ch_per_region = num_ch_per_region;
            STOPmrv_metadata.num_regions_complete = iRegion;

%             save(vecDiffmat_saveName, 're_mean_mrv', 'im_mean_mrv', 'STOPmrv_metadata');
            
        end    % for iRegion...
        
    end    % for iSession...
    
end
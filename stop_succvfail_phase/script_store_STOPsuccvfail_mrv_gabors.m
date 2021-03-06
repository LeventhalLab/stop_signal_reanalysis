% script_store_STOPsuccvfail_mrv_gabors

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

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
    
    for iSession = 1 : numSessions                 % CHANGE THIS BACK TO 1 TO REDO ALL OF THE CALCULATIONS
        
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
            stopPhase_name = ['stopPhases_' sessionChannels{iCh}.name '_gabor.mat'];
            stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
            if ~exist(stopPhase_name, 'file'); continue; end
            
            load(stopPhase_name);
            break;
        end
        if ~exist(stopPhase_name, 'file'); continue; end
        STOPmrv_metadata.Fs = STOPmetadata.Fs;
        STOPmrv_metadata.f = STOPmetadata.f;
        STOPmrv_metadata.eventList = STOPmetadata.eventList;
        STOPmrv_metadata.chName = STOPmetadata.chName;
        STOPmrv_metadata.twin = STOPmetadata.twin;
        
        numEventTypes = length(STOPmetadata.eventList);
        numFreqs  = length(STOPmetadata.f);
        numSamps  = size(correctSTOP_phase, 2);

%         mrvDiff = zeros(numEventTypes, numFreqs, numSamps);
%         mrv = zeros(2, numEventTypes, numFreqs, numSamps);
        mrv = zeros(2, numEventTypes, numSamps, numFreqs);
        t = linspace(STOPmrv_metadata.twin(1),STOPmrv_metadata.twin(2),numSamps);
        f = STOPmetadata.f;
        
        STOPmrv_metadata.t = t;
        STOPmrv_metadata.regionList = regionList;
        STOPmrv_metadata.chNames = cell(1, numRegions);
        
        vecDiffmat_saveName = ['vecDiff_stopSuccvFail_' sessionList{iSession} '_gabors.mat'];
        vecDiffmat_saveName = fullfile(session_stopPhaseDir, vecDiffmat_saveName);
        if exist(vecDiffmat_saveName,'file')
            load(vecDiffmat_saveName);
            if STOPmrv_metadata.num_regions_complete == numRegions; continue; end
            num_ch_per_region = STOPmrv_metadata.num_ch_per_region;
        else
%             mrvDiff_by_region = cell(1, numRegions);
            mrv_by_region = cell(1, numRegions);
%             mean_mrvDiff = zeros(numRegions,numEventTypes,numFreqs,numSamps);    % left over from Hilbert version
            mean_mrvDiff = zeros(numRegions, numEventTypes, numSamps, numFreqs);
            mean_mrv = zeros(numRegions, 2, numEventTypes, numSamps, numFreqs);   % 1 for correct, 2 for failed STOP
            num_ch_per_region = zeros(1, numRegions);
            STOPmrv_metadata.num_regions_complete = 0;
        end

        for iRegion = STOPmrv_metadata.num_regions_complete + 1 : numRegions
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numRegionChannels = length(regionChannels);
            
%             mrvDiff_by_region{iRegion} = zeros(numRegionChannels, numEventTypes, numFreqs, numSamps);
            mrv_by_region{iRegion} = zeros(numRegionChannels, 2, numEventTypes, numSamps, numFreqs);
            
            num_ch_per_region(iRegion) = numRegionChannels;
            STOPmrv_metadata.chNames{iRegion} = cell(1, numRegionChannels);

            for iCh = 1 : numRegionChannels
%                 iCh
            
                ch = regionChannels{iCh};
                STOPmrv_metadata.chNames{iRegion}{iCh} = ch.name;
                mrv_mat_saveName = ['mrv_stopSuccvFail_' ch.name '_gabor.mat'];
                mrv_mat_saveName = fullfile(session_stopPhaseDir, mrv_mat_saveName);
                if exist(mrv_mat_saveName, 'file'); continue; end
            
                stopPhase_name = ['stopPhases_' ch.name '_gabor.mat'];
                stopPhase_name = fullfile(session_stopPhaseDir, stopPhase_name);
                if ~exist(stopPhase_name, 'file'); continue; end

                surr_mrv_mat_saveName = ['mrv_stopSuccvFail_surrogates_' ch.name '_gabor.mat'];
                surr_mrv_mat_saveName = fullfile(session_stopPhaseDir, surr_mrv_mat_saveName);
                if ~exist(surr_mrv_mat_saveName, 'file'); continue; end
                
                load(stopPhase_name);
                load(surr_mrv_mat_saveName);
                
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

            save(vecDiffmat_saveName, 're_mean_mrv', 'im_mean_mrv', 'STOPmrv_metadata');
            
        end    % for iRegion...
        
    end    % for iSession...
    
end
% script_plot_STOPsuccvfail_mrv_z_across_sessions_c

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

regions_per_page = 5;
colorLim = [0 1];
vecDiff_colorLim = [0 1];
z_colorLim = [-3 3];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

for i_chDB = 3 : 3%length(chDB_list)
    
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
    
    cp = initChanParams();
    cp.task = 3;    % only stop-signal sessions
    chList = extractChannels( cp, channels );
    STOPchannels = channels( chList );
        
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
    STOPchannels = excludeChannels(cp, STOPchannels);

    sessionList = getSessionsfromChannelDB( STOPchannels );
    numSessions = length(sessionList);
    
    allRegionList = getRegionsfromChannelDB( STOPchannels );
    num_allRegions = length(allRegionList);
    
    for i_startSession = 1 : numSessions
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{i_startSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{i_startSession} '...'])
            continue
        end
        z_DiffmatName = ['region_z_stopSuccvFail_' sessionList{i_startSession} '.mat'];
        z_DiffmatName = fullfile(session_stopPhaseDir, z_DiffmatName);
        if ~exist(z_DiffmatName,'file'); continue; end
        load(z_DiffmatName);
        
        f = z_STOPmrv_metadata.f; t = z_STOPmrv_metadata.t;
        numEventTypes = length(z_STOPmrv_metadata.eventList);
        numFreqs  = length(f);
        numSamps  = length(t);
        
        break;
    end
        
    region_z = zeros(num_allRegions, numEventTypes, numFreqs, numSamps);
    numRegionSessions = zeros(1, num_allRegions);
    
    for iSession = i_startSession : numSessions
        disp(sprintf('Session %d out of %d', iSession, numSessions));
        if iSession > i_startSession    % already loaded vecDiffmatName above
            session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
            if ~exist(session_stopPhaseDir, 'dir')
                disp([session_stopPhaseDir ' not found. Skipping ' sessionList{iSession} '...'])
                continue
            end
            z_DiffmatName = ['region_z_stopSuccvFail_' sessionList{iSession} '.mat'];
            z_DiffmatName = fullfile(session_stopPhaseDir, z_DiffmatName);
            if ~exist(z_DiffmatName,'file'); continue; end
            load(z_DiffmatName);
        end
        
        STOPmrv_z_acrossSessions_metadata(iSession) = z_STOPmrv_metadata;
        numSessionRegions = length(z_STOPmrv_metadata.regionList);
        for iSessionRegion = 1 : numSessionRegions
            
            allRegionIdx = find(strcmpi(z_STOPmrv_metadata.regionList{iSessionRegion}, ...
                                        allRegionList));
%             mean_z = re_mean_z + 1i*im_mean_z;                        
            region_z(allRegionIdx,:,:,:) = region_z(allRegionIdx,:,:,:) + ...
                                           mean_z(iSessionRegion,:,:,:);
            numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
            
        end

    end
    
    figProps.n = numEventTypes;
    figProps.colSpacing = ones(1, figProps.n) * 0.5;
    figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                 2 * sideMargins - ...
                                                 sum(figProps.colSpacing)) / figProps.n;

    zphase_acrossSessions_fig_saveName = [implantID '_zphase_stopSuccvFail.pdf'];
    zphase_acrossSessions_fig_saveName = fullfile(subject_stopPhaseDir, zphase_acrossSessions_fig_saveName);
    
    vecDiff_z_acrossSessions_mat_saveName = [implantID '_vecDiff_z_stopSuccvFail.mat'];
    vecDiff_z_acrossSessions_mat_saveName = fullfile(subject_stopPhaseDir, vecDiff_z_acrossSessions_mat_saveName);
        
    numRegionPlots = 0;
    numPages = 0;
    for iRegion = 1 : num_allRegions
%         cp = initChanParams();
%         cp.locationName = allRegionList{iRegion};            
%         chList = extractChannels(cp, sessionChannels);
%         if isempty(chList); continue; end
%         regionChannels = sessionChannels(chList);
%         numRegionChannels = length(regionChannels);

        region_z(iRegion,:,:,:) = region_z(iRegion,:,:,:) / numRegionSessions(iRegion);

%         mrvDiff = squeeze(region_mrv(iRegion,1,:,:,:) - region_mrv(iRegion,2,:,:,:));
                
        numRegionPlots = numRegionPlots + 1;
        rowNum = rem(numRegionPlots, regions_per_page);
        if rowNum == 1
            [h_vecDiff_fig, h_vecDiff_axes] = createFigPanels5(figProps);
            page_regionList = allRegionList{iRegion};   %ch.name;
            page_numSessionList = num2str(numRegionSessions(iRegion));
            numPages = numPages + 1;
        else
            page_regionList = [page_regionList ', ' allRegionList{iRegion}];
            page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
        end
        if rowNum == 0; rowNum = regions_per_page; end

        for iEventType = 1 : numEventTypes
            axes(h_vecDiff_axes(rowNum, iEventType));

            z_toPlot = squeeze(region_z(iRegion,iEventType, :, :));
            imagesc(t,f,z_toPlot);
            set(gca,'ydir','normal','clim', z_colorLim);
            if rowNum == 1
                title(z_STOPmrv_metadata.eventList{iEventType});
            end
            if rowNum < regions_per_page
                set(gca,'xticklabel',[]);
            end
            if iEventType > 1
                set(gca,'yticklabel',[]);
            end
        end

        if (rowNum == regions_per_page || iRegion == num_allRegions)

            h_figAxes = createFigAxes(h_vecDiff_fig);
            axes(h_figAxes);

            textStr{1} = ['difference in mean resultant vectors (z-scores), averaged across sessions'];
            textStr{2} = page_regionList;
            textStr{3} = ['number of sessions per region: ' page_numSessionList];
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            if numPages == 1
                export_fig(zphase_acrossSessions_fig_saveName,'-pdf','-q101','-painters');
            else
                export_fig(zphase_acrossSessions_fig_saveName,'-pdf','-q101','-painters','-append');
            end
            close(h_vecDiff_fig);

        end

    end    % for iRegion...
    
    save(vecDiff_z_acrossSessions_mat_saveName, 'region_z', 'allRegionList', 'numRegionSessions', 'STOPmrv_z_acrossSessions_metadata');
        
end    % for i_chDB...


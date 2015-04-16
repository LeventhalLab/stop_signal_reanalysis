% script_plot_STOPsuccvfail_mrv_surrDiff

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

channels_per_page = 5;
colorLim = [0 1];
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

for i_chDB = 2 : 2%length(chDB_list)
    
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
    
    for iSession = 1 : numSessions
        
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
        
        tic
        for iCh = 1 : numSessionChannels
            mrv_matName = ['mrv_stopSuccvFail_' sessionChannels{iCh}.name '.mat'];
            mrv_matName = fullfile(session_stopPhaseDir, mrv_matName);
            if ~exist(mrv_matName, 'file'); continue; end
            
            load(mrv_matName);
            break;
        end
        toc

        f = STOPmrv_metadata.f; t = STOPmrv_metadata.t;
        numEventTypes = length(STOPmrv_metadata.eventList);
        numFreqs  = length(f);
        numSamps  = length(t);
        
        z_STOPmrv_metadata.Fs = STOPmrv_metadata.Fs;
        z_STOPmrv_metadata.freqBands = STOPmrv_metadata.freqBands;
        z_STOPmrv_metadata.eventList = STOPmrv_metadata.eventList;
        z_STOPmrv_metadata.twin = STOPmrv_metadata.twin;
        z_STOPmrv_metadata.regionList = regionList;
        z_STOPmrv_metadata.f = f; z_STOPmrv_metadata.t = t;
        
        mean_z = zeros(numRegions, numEventTypes, numFreqs, numSamps);
        ch_per_region = zeros(1, numRegions);
            
        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
                                                 
        surr_vecDiff_fig_saveName = ['z_stopSuccvFail_' sessionList{iSession} '.pdf'];
        surr_vecDiff_fig_saveName = fullfile(session_stopPhaseDir, surr_vecDiff_fig_saveName);
        
%         vecDiffmatName = ['vecDiff_stopSuccvFail_' sessionList{iSession} '.mat'];
%         vecDiffmatName = fullfile(session_stopPhaseDir, vecDiffmatName);
        
        z_DiffmatName = ['region_z_stopSuccvFail_' sessionList{iSession} '.mat'];
        z_DiffmatName = fullfile(session_stopPhaseDir, z_DiffmatName);
        
%         if ~exist(vecDiffmatName,'file'); continue; end

        numChPlots = 0;
        numPages = 0;

        for iRegion = 1 : numRegions
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numRegionChannels = length(regionChannels);
            ch_per_region(iRegion) = numRegionChannels;
            
            for iCh = 1 : numRegionChannels
                iCh

%                 STOPmrv_metadata.eventsComplete = 0;
%                 STOPmrv_metadata_loaded = false;
            
                ch = regionChannels{iCh};
                mrv_matName = ['mrv_stopSuccvFail_' ch.name '.mat'];
                mrv_matName = fullfile(session_stopPhaseDir, mrv_matName);
                if ~exist(mrv_matName, 'file'); continue; end
                load(mrv_matName);
                mrv = re_mrv + 1i*im_mrv;
                
                mrv_surrogate_name = ['mrv_stopSuccvFail_surrogates_' ch.name '.mat'];
                mrv_surrogate_name = fullfile(session_stopPhaseDir, mrv_surrogate_name);
                if ~exist(mrv_matName, 'file'); continue; end
                load(mrv_surrogate_name);
                
                mrvDiff = squeeze(mrv(1,:,:,:) - mrv(2,:,:,:));
                z = (abs(mrvDiff) - squeeze(mean_surr_diff)) ./ squeeze(std_surr_diff);
                mean_z(iRegion, :, :, :) = squeeze(mean_z(iRegion, :, :, :)) + z;
                
                numChPlots = numChPlots + 1;
                rowNum = rem(numChPlots, channels_per_page);
                if rowNum == 1
                    [h_vecDiff_fig, h_vecDiff_axes] = createFigPanels5(figProps);
                    page_chList = STOPmrv_metadata.chNames{iRegion}{iCh};   %ch.name;
                    page_locList = STOPmrv_metadata.regionList{iRegion};    %ch.location.subclass;
                    numPages = numPages + 1;
                else
                    page_chList = [page_chList ', ' ch.name];
                    page_locList = [page_locList ', ' ch.location.subclass];
                end
                if rowNum == 0; rowNum = channels_per_page; end

                for iEventType = 1 : numEventTypes
                    axes(h_vecDiff_axes(rowNum, iEventType));
                    
                    z_toPlot = squeeze(z(iEventType, :, :));
                    imagesc(t,f,z_toPlot);
                    set(gca,'ydir','normal','clim', z_colorLim);
                    if rowNum == 1
                        title(STOPmrv_metadata.eventList{iEventType});
                    end
                    if rowNum < channels_per_page
                        set(gca,'xticklabel',[]);
                    end
                    if iEventType > 1
                        set(gca,'yticklabel',[]);
                    end
                end

                if (rowNum == channels_per_page || numChPlots == numSessionChannels)

                    h_figAxes = createFigAxes(h_vecDiff_fig);
                    axes(h_figAxes);

                    textStr{1} = ['z-score, difference in mean resultant vectors'];
                    textStr{2} = page_chList;
                    textStr{3} = page_locList;
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                    if numPages == 1
                        export_fig(surr_vecDiff_fig_saveName,'-pdf','-q101','-painters');
                    else
                        export_fig(surr_vecDiff_fig_saveName,'-pdf','-q101','-painters','-append');
                    end
                    close(h_vecDiff_fig);

                end
                    
            end    % for iCh...
            
            mean_z(iRegion, :, :, :) = mean_z(iRegion, :, :, :) / ch_per_region(iRegion);      

        end    % for iRegion...
        
        save(z_DiffmatName, 'mean_z', 'ch_per_region', 'z_STOPmrv_metadata');
        
        % create a new page of figures to show the region means for
        % each session
%         load(vecDiffmatName);
%         mean_mrv = re_mean_mrv + 1i*im_mean_mrv;
%         mean_mrvDiff = squeeze(mean_mrv(:, 1, :, :, :) - mean_mrv(:, 2, :, :, :));
        for iRegion = 1 : numRegions                                                        % WORKING HERE...
            rowNum = rem(iRegion, channels_per_page);
            if rowNum == 1
                [h_vecDiff_fig, h_vecDiff_axes] = createFigPanels5(figProps);
                page_regionList = regionList{iRegion};
                page_numChList  = ['Number of Channels per Region: ' num2str(ch_per_region(iRegion))];
            else
                page_regionList = [page_regionList ', ' regionList{iRegion}];
                page_numChList  = [page_numChList ', ' num2str(ch_per_region(iRegion))];
            end
            if rowNum == 0; rowNum = channels_per_page; end

            for iEventType = 1 : numEventTypes
                axes(h_vecDiff_axes(rowNum, iEventType));

                toPlot = (squeeze(mean_z(iRegion, iEventType, :, :)));
                imagesc(t,f,toPlot);
                set(gca,'ydir','normal','clim',z_colorLim);
                if rowNum == 1
                    title(STOPmrv_metadata.eventList{iEventType});
                end
                if rowNum < channels_per_page
                    set(gca,'xticklabel',[]);
                end
                if iEventType > 1
                    set(gca,'yticklabel',[]);
                end
            end
            if (rowNum == channels_per_page || iRegion == numRegions)
                h_figAxes = createFigAxes(h_vecDiff_fig);
                axes(h_figAxes);

                textStr{1} = ['mean difference in mean resultant vectors (z-score) by region, stop-success vs stop-fail'];
                textStr{2} = page_regionList;
                textStr{3} = page_numChList;
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                export_fig(surr_vecDiff_fig_saveName,'-pdf','-q101','-painters','-append');

                close(h_vecDiff_fig);
            end

        end    % for iRegion...


    end    % for iSession...
    
end
% script_plotphaseRTcorrelations

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

makePlots = true;

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
channels_per_page = 5;
colorLim = [0 1];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
                                                  
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
    
%     regionList = getRegionsfromChannelDB( channels );
%     numRegions = length(regionList);
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        regionList = getRegionsfromChannelDB(sessionChannels);
        numRegions = length(regionList);
        
        totalSessionChannels = length(sessionChannels);
        
        phaseRThist_sessionDir = fullfile(subject_phaseRThistdir, sessionList{iSession});
        if ~exist(phaseRThist_sessionDir, 'dir')
            continue;
        end
        
        fig_saveName = ['phase_RT_circCorrPlots_' sessionList{iSession} '.pdf'];
        fig_saveName = fullfile(phaseRThist_sessionDir, fig_saveName);
%         if exist(saveName,'file'); continue; end

        region_mean_saveName = ['phase_RTbyRegion_' sessionList{iSession} '.mat'];
        region_mean_saveName = fullfile(phaseRThist_sessionDir, region_mean_saveName);
        
        
        if exist(region_mean_saveName, 'file'); continue; end  % comment this line out to make plots whether means have been saved or not
        
        
        
        corr_fileName = ['phase_RT_circCorr_' sessionChannels{1}.name '.mat'];
        corr_fileName = fullfile(phaseRThist_sessionDir, corr_fileName);
        if ~exist(corr_fileName, 'file'); continue; end
        load( corr_fileName );
        
        numEventTypes = length(phaseRThist_metadata.eventList);
        Fs = phaseRThist_metadata.Fs;
        numSamps = size(pcc, 2);
        tmin = phaseRThist_metadata.twin(1);
        
        t = linspace(tmin + 1/Fs, tmin + numSamps/Fs, numSamps);
        f = phaseRThist_metadata.freqs;
        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
        numPages = 0;
        pcc_by_region = cell(1, numRegions);
        mean_pcc = zeros(numRegions, numEventTypes, length(t), length(f));
        num_ch_per_region = zeros(1, numRegions);
        numChPlots = 0;
        for iRegion = 1 : numRegions
            
            cp = initChanParams();
            cp.locationName = regionList{iRegion};                  % DEBUGGING HERE
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numCh = length(regionChannels);
            
            pcc_by_region{iRegion} = zeros(numCh, numEventTypes, length(t), length(f));
            num_ch_per_region(iRegion) = numCh;
            
            for iCh = 1 : numCh
    %             iCh
                ch = regionChannels{iCh};
                corr_fileName = ['phase_RT_circCorr_' ch.name '.mat'];
                corr_fileName = fullfile(phaseRThist_sessionDir, corr_fileName);
                if ~exist(corr_fileName, 'file'); continue; end

                load( corr_fileName );

                numChPlots = numChPlots + 1;
                rowNum = rem(numChPlots, channels_per_page);
                if rowNum == 1
                    if makePlots
                        [h_fig, h_axes] = createFigPanels5(figProps);
                    end
                    page_chList = ch.name;
                    page_locList = ch.location.subclass;
                    numPages = numPages + 1;
                else
                    page_chList = [page_chList ', ' ch.name];
                    page_locList = [page_locList ', ' ch.location.subclass];
                end
                if rowNum == 0; rowNum = channels_per_page; end

                pcc_by_region{iRegion}(iCh, :, :, :) = pcc;
                
                for iEventType = 1 : numEventTypes
                    toPlot = squeeze(pcc(iEventType, :, :))';

                    if makePlots
                        axes(h_axes(rowNum, iEventType));
                        imagesc(t, f, toPlot);
                        set(gca,'ydir','normal','clim',colorLim);
                        if rowNum == 1
                            title(phaseRThist_metadata.eventList{iEventType});
                        end
                        if rowNum < channels_per_page
                            set(gca,'xticklabel',[]);
                        end
                        if iEventType > 1
                            set(gca,'yticklabel',[]);
                        end
                    end
                end

                if (rowNum == channels_per_page || numChPlots == totalSessionChannels) && makePlots
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);

                    textStr{1} = page_chList;
                    textStr{2} = page_locList;
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                    if numPages == 1
                        export_fig(fig_saveName,'-pdf','-q101','-painters');
                    else
                        export_fig(fig_saveName,'-pdf','-q101','-painters','-append');
                    end
                    close(h_fig);
                end

            end    % for iCh...
            % create averages within individual brain regions for each session
            mean_pcc(iRegion, : ,:, :) = mean(pcc_by_region{iRegion}, 1);
            
        end    % for iRegion...
        region_phaseRTcorr_metadata = phaseRThist_metadata;
        region_phaseRTcorr_metadata.regionList = regionList;
        region_phaseRTcorr_metadata.num_ch_per_region = num_ch_per_region;
        save(region_mean_saveName, 'mean_pcc', 'region_phaseRTcorr_metadata');
        
        % create a new page of figures to show the region means for
        % each session
        if makePlots
            for iRegion = 1 : numRegions
                rowNum = rem(iRegion, channels_per_page);
                if rowNum == 1
                    [h_fig, h_axes] = createFigPanels5(figProps);
                    page_regionList = regionList{iRegion};
                    page_numChList  = ['Number of Channels per Region: ' num2str(num_ch_per_region(iRegion))];
                else
                    page_regionList = [page_regionList ', ' regionList{iRegion}];
                    page_numChList  = [page_numChList ', ' num2str(num_ch_per_region(iRegion))];
                end
                if rowNum == 0; rowNum = channels_per_page; end
                
                for iEventType = 1 : numEventTypes
                    toPlot = squeeze(mean_pcc(iRegion, iEventType, :, :))';
                    axes(h_axes(rowNum, iEventType));
                    imagesc(t, f, toPlot);
                    set(gca,'ydir','normal','clim',colorLim);
                    if rowNum == 1
                        title(phaseRThist_metadata.eventList{iEventType});
                    end
                    if rowNum < channels_per_page
                        set(gca,'xticklabel',[]);
                    end
                    if iEventType > 1
                        set(gca,'yticklabel',[]);
                    end
                end

                if (rowNum == channels_per_page || iRegion == numRegions)
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);

                    textStr{1} = 'Average phase-RT Correlations within Regions';
                    textStr{2} = page_regionList;
                    textStr{3} = page_numChList;
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                    export_fig(fig_saveName,'-pdf','-q101','-painters','-append');
                    close(h_fig);
                end
            end    % for iRegion...
        end    % if makePlots...
                    
    end    % end for iSession...
    
end    % for i_chDB...
                
            
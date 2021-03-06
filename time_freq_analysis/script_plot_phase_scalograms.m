% script to plot phase scalograms

chDB_directory    = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
scalogramDir = '/Volumes/Tbolt_02/stop-signal reanalysis/trial_scalograms';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};
numTrialTypes = length(trialTypeList);

regions_per_page = 5;
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;
figProps.m = regions_per_page;    % number of rows
sideMargins = 2; botMargin = 2.54;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

colLim = [0 0.5];
desired_freq_ticks = [2,8,20,50,80,250];

cmap = 'jet';

for i_chDB = 1:4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    
    subject_scalogramDir = fullfile(scalogramDir, [implantID '_ps']);
    if ~exist(subject_scalogramDir,'dir'); continue; end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iTrialType = 1 : numTrialTypes
        trialType = trialTypeList{iTrialType};
        
        for iCh = 1 : length(channels)
            ch = channels{iCh};
            session_scalogramDir = fullfile(subject_scalogramDir,[ch.session '_scalograms']);
            test_ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
            if ~exist(test_ch_scalogramName,'file');continue;end
            break;
        end
        load(test_ch_scalogramName);
        numEvents = length(scalogram_metadata.eventList);
        figProps.n = numEvents;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
        t = scalogram_metadata.t; f = scalogram_metadata.f;
        f_idx = 1:length(f);
        freqTick_idx = zeros(1, length(desired_freq_ticks));
        freqTick_label = zeros(1, length(desired_freq_ticks));
        for i_freqTick = 1 : length(desired_freq_ticks)
            freqTick_idx(i_freqTick) = find(abs(f - desired_freq_ticks(i_freqTick)) == ...
                                            min(abs(f - desired_freq_ticks(i_freqTick))));
            freqTick_label(i_freqTick) = round(f(freqTick_idx(i_freqTick)));
        end
                                        

        numSamps  = length(t); numFreqs = length(f);
        mean_sessionRegionMRL = NaN(numSessions, numRegions, numEvents, numSamps, numFreqs);
            
        numPages = 0;
        for iSession = 1 : numSessions
            session_scalogramDir = fullfile(subject_scalogramDir,[sessionList{iSession} '_scalograms']);
            if ~exist(session_scalogramDir,'file'); continue; end

            saveName = [sessionList{iSession} '_' trialType '_mrl'];
            mean_sessionMRLname = fullfile(session_scalogramDir,[sessionList{iSession} '_' trialType '_meanRegionMRL.mat']);
            
            fprintf('session %s, %d of %d\n', ...
                sessionList{iSession}, iSession, numSessions)

            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.locationSubClass = ROI_list;

            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
            if isempty(sessionChannels);continue;end
            
            numSessionChannels = length(sessionChannels);

            for iRegion = 1 : numRegions
                cp = initChanParams();
                cp.locationSubClass = ROI_list{iRegion};
                region_chList = extractChannels( cp, sessionChannels );
                sessionRegionChannels = sessionChannels(region_chList);
                if isempty(sessionRegionChannels);continue;end
                numSessionRegionChannels = length(region_chList);
                
                mean_chMRL = NaN(numSessionRegionChannels, numEvents, numSamps, numFreqs);
                
                numChPlots = 0;
                figProps.m = regions_per_page;
                figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                              figProps.topMargin - botMargin - ...
                                                              sum(figProps.rowSpacing)) / figProps.m;
                for iCh = 1 : numSessionRegionChannels
                    ch = sessionRegionChannels{iCh};
                    ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
                    mean_chMRLname = fullfile(session_scalogramDir,[ch.name '_' trialType '_meanMRL.mat']);
                    if numSessionRegionChannels < figProps.m
                        figProps.m = numSessionRegionChannels;
                        figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                                      figProps.topMargin - botMargin - ...
                                                                      sum(figProps.rowSpacing)) / figProps.m;
                    end
                    if ~exist(ch_scalogramName,'file');continue;end
                    if ~strcmpi(ch_scalogramName, test_ch_scalogramName)
                        load(ch_scalogramName);
                    end
                    
                    t_ticks = [scalogram_metadata.twin(1),0,scalogram_metadata.twin(2)];
                    
                    numChPlots = numChPlots + 1;
                    plotRow = rem(numChPlots,figProps.m);
                    if plotRow == 0
                        plotRow = figProps.m;
                    end
                    if plotRow == 1
                        [h_fig,h_axes] = createFigPanels5(figProps);
                        page_chList = ch.name;
                        numPages = numPages + 1;
                    else
                        page_chList = [page_chList ', ' ch.name];
                    end

                    Wphase = angle(W);
                    mrl = abs(squeeze(mean(exp(1i * Wphase),3)));
                    
                    if length(scalogram_metadata.eventList) == 1
                        mean_chMRL(iCh,1,:,:) = mrl;
                    else
                        mean_chMRL(iCh,:,:,:) = mrl;
                    end
                    
                    if ~exist(mean_chMRLname, 'file')
                        save(mean_chMRLname, 'mean_chMRL', 'scalogram_metadata');
                    end
                    
                    for iEvent = 1 : length(scalogram_metadata.eventList)
                        axes(h_axes(plotRow, iEvent));
                        
                        toPlot = squeeze(mean_chMRL(iCh,iEvent,:,:))';
                        h_pcolor = pcolor(scalogram_metadata.t, ...
                                f, ...
                                toPlot);
                        h_pcolor.EdgeColor = 'none';
                        colormap(cmap);
                        set(gca,'yscale','log');
                        set(gca,'ydir','normal',...
                                'xtick',t_ticks,...
                                'clim',colLim,...
                                'ytick',f(freqTick_idx));%f_idx(1:16:length(f)));
                        
                        if plotRow == 1
                            title(scalogram_metadata.eventList{iEvent});
                        end
                        if iEvent > 1
                            set(gca,'yticklabel',[]);
                        else
                            set(gca,'yticklabel',freqTick_label);%round(f(1:16:length(f))));
                        end
                        
                        if plotRow < figProps.m
                            set(gca,'xticklabel',[]);
                        end
                    end    % for iEvent...
                    
                    if plotRow == figProps.m || iCh == numSessionRegionChannels
                        h_figAxes = createFigAxes(h_fig);
                        axes(h_figAxes);
                        
                        cur_PDFname = sprintf('%s_%02d.pdf', saveName, numPages);
                        cur_figName = sprintf('%s_%02d.fig', saveName, numPages);
                        cur_PDFname = fullfile(session_scalogramDir, cur_PDFname);
                        cur_figName = fullfile(session_scalogramDir, cur_figName);
                        
                        textStr = cell(1,3);
                        textStr{1} = [implantID ', ' trialType ', MRL, single channels in ' ROI_list{iRegion}];
                        textStr{2} = page_chList;
                        textStr{3} = ['color limits: ' num2str(colLim(1)) ', ' num2str(colLim(2))];
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                        print(cur_PDFname, '-dpdf');
                        savefig(h_fig,cur_figName,'compact');
%                         if numPages == 1
%                             export_fig(PDFname, '-pdf', '-q101', '-painters','-nocrop');
%                         else
%                             export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
%                         end
                        close(h_fig);
                    end
                end    % for iCh...
                mean_sessionRegionMRL(iSession, iRegion, :, :, :) = squeeze(mean(mean_chMRL,1));
                
            end    % for iRegion...
            % now make a page of plots with the average for each region for
            % the current session
            figProps.m = regions_per_page;
            figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                          figProps.topMargin - botMargin - ...
                                                          sum(figProps.rowSpacing)) / figProps.m;
            for iRegion = 1 : numRegions
                if numRegions < figProps.m
                    figProps.m = numRegions;
                    figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                                  figProps.topMargin - botMargin - ...
                                                                  sum(figProps.rowSpacing)) / figProps.m;
                end
                
                plotRow = rem(iRegion,figProps.m);
                if plotRow == 0
                    plotRow = figProps.m;
                end
                if plotRow == 1
                    [h_fig,h_axes] = createFigPanels5(figProps);
                    page_regionList = ROI_list{iRegion};
                    numPages = numPages + 1;
                else
                    page_regionList = [page_regionList ', ' ROI_list{iRegion}];
                end
                
                for iEvent = 1 : numEvents
                    axes(h_axes(plotRow, iEvent));
                    
                    toPlot = squeeze(mean_sessionRegionMRL(iSession,iRegion,iEvent,:,:))';
                    h_pcolor = pcolor(scalogram_metadata.t, ...
                            f, ...
                            toPlot);
                    h_pcolor.EdgeColor = 'none';
                    colormap(cmap);
                    set(gca,'yscale','log');
                    set(gca,'ydir','normal',...
                            'xtick',t_ticks,...
                            'clim',colLim,...
                            'ytick',f(freqTick_idx));%f_idx(1:16:length(f)));
                        
                    if plotRow == 1
                        title(scalogram_metadata.eventList{iEvent});
                    end
                    if iEvent > 1
                        set(gca,'yticklabel',[]);
                    else
                        set(gca,'yticklabel',freqTick_label);%round(f(1:16:length(f))));
                    end

                    if plotRow < figProps.m
                        set(gca,'xticklabel',[]);
                    end
                end    % for iEvent...
                
                if plotRow == figProps.m || iRegion == numRegions
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);
                    
                    cur_PDFname = sprintf('%s_region_means_%02d.pdf', saveName, numPages);
                    cur_figName = sprintf('%s_region_means_%02d.fig', saveName, numPages);
                    cur_PDFname = fullfile(session_scalogramDir, cur_PDFname);
                    cur_figName = fullfile(session_scalogramDir, cur_figName);
                        
                    textStr = cell(1,4);
                    textStr{1} = [implantID ', ' trialType ', MRL, averaged across channels for each region'];
                    textStr{2} = sessionList{iSession};
                    textStr{3} = page_regionList;
                    textStr{4} = ['color limits: ' num2str(colLim(1)) ', ' num2str(colLim(2))];
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                    print(cur_PDFname, '-dpdf');
                    savefig(h_fig,cur_figName,'compact');

                    close(h_fig);
                end
                
            end    % for iRegion...
                
        end    % for iSession...
        
        % now, average across all session-regions
        saveName = [implantID '_' trialType '_mrl_scalograms'];
        mean_subjectPowerName = fullfile(subject_scalogramDir,[implantID '_' trialType '_meanScalograms.mat']);
        figProps.m = regions_per_page;
        figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                      figProps.topMargin - botMargin - ...
                                                      sum(figProps.rowSpacing)) / figProps.m;
        for iRegion = 1 : numRegions
            if numRegions < figProps.m
                figProps.m = numRegions;
                figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                              figProps.topMargin - botMargin - ...
                                                              sum(figProps.rowSpacing)) / figProps.m;
            end
            
            meanRegionMRL = squeeze(nanmean(mean_sessionRegionMRL(:,iRegion,:,:,:),1));

            plotRow = rem(iRegion,figProps.m);
            if plotRow == 0
                plotRow = figProps.m;
            end
            if plotRow == 1
                [h_fig,h_axes] = createFigPanels5(figProps);
                page_regionList = ROI_list{iRegion};
                numPages = numPages + 1;
            else
                page_regionList = [page_regionList ', ' ROI_list{iRegion}];
            end
            
            for iEvent = 1 : numEvents
                axes(h_axes(plotRow, iEvent));
                
                if numEvents == 1
                    toPlot = meanRegionMRL';
                else
                    toPlot = squeeze(meanRegionMRL(iEvent,:,:))';
                end
                h_pcolor = pcolor(scalogram_metadata.t, ...
                        f, ...
                        toPlot);
                h_pcolor.EdgeColor = 'none';
                colormap(cmap);
                set(gca,'yscale','log');
                set(gca,'ydir','normal',...
                        'xtick',t_ticks,...
                        'clim',colLim,...
                        'ytick',f(freqTick_idx));%round(f(1:16:length(f))));

                if plotRow == 1
                    title(scalogram_metadata.eventList{iEvent});
                end
                if iEvent > 1
                    set(gca,'yticklabel',[]);
                else
                    set(gca,'yticklabel',freqTick_label);%f_idx(1:16:length(f)));
                end

                if plotRow < figProps.m
                    set(gca,'xticklabel',[]);
                end
            end    % for iEvent...

            if plotRow == figProps.m || iRegion == numRegions
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);
                
                cur_PDFname = sprintf('%s_%02d.pdf', saveName, numPages);
                cur_figName = sprintf('%s_%02d.fig', saveName, numPages);
                cur_PDFname = fullfile(subject_scalogramDir, cur_PDFname);
                cur_figName = fullfile(subject_scalogramDir, cur_figName);
                        
                textStr = cell(1,3);
                textStr{1} = [implantID ', ' trialType ', MRL, averaged across channels and sessions for each region'];
                textStr{2} = page_regionList;
                textStr{3} = ['color limits: ' num2str(colLim(1)) ', ' num2str(colLim(2))];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                print(cur_PDFname, '-dpdf');
                savefig(h_fig,cur_figName,'compact');
                close(h_fig);
            end
            
        end    % for iRegion...

        subject_region_metadata.sessionList = sessionList;
        subject_region_metadata.Fs = scalogram_metadata.Fs;
        subject_region_metadata.trialType = scalogram_metadata.trialType;
        subject_region_metadata.t = scalogram_metadata.t;
        subject_region_metadata.f = scalogram_metadata.f;
        subject_region_metadata.twin = scalogram_metadata.twin;
        subject_region_metadata.eventList = scalogram_metadata.eventList;
        subject_region_metadata.regionList = ROI_list;
        % may need more metadata, can add them later
        save(mean_subjectPowerName, 'mean_sessionRegionMRL', 'subject_region_metadata');
        
    end    % for iTrialType...
    
end    % for i_chDB...



            




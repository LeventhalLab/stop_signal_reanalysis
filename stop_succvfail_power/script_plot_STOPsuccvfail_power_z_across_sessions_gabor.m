% script_plot_STOPsuccvfail_mrv_z_across_sessions_gabor

% script to go through all STOP sessions and see if there are consistent
% power differences between STOP-success and STOP-failure trials

% chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
% power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power_gabors';

chDB_directory    = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
phase_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

y_zoom_lim = [1 100];   % for more focused plot
y_full_plot_lim = [0.1 502];

regions_per_page = 5;
colorLim = [0 1];
power_colorLim = [0 1];
z_colorLim = [-2.5 2.5];
% desired_freq_ticks = [4,8,16,32,64,128,256,500];
desired_freq_ticks = [4,16,64,128,500];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = regions_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

for i_chDB = 1 : 4%length(chDB_list)
    
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
    cp.task = 3;    % only stop-signal sessions
    chList = extractChannels( cp, channels );
    STOPchannels = channels( chList );
        
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
    STOPchannels = excludeChannels(cp, STOPchannels);
    
    cp = initChanParams();
    cp.locationSubClass = {'e2', 'e3', 'e02','e03'};
    STOPchannels = excludeChannels(cp, STOPchannels);
    if isempty(STOPchannels);continue;end

    sessionList = getSessionsfromChannelDB( STOPchannels );
    numSessions = length(sessionList);
    
    allRegionList = getRegionsfromChannelDB( STOPchannels );
    num_allRegions = length(allRegionList);
    
%     for i_startSession = 1 : numSessions
%         
%         session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{i_startSession});
%         if ~exist(session_stopPhaseDir, 'dir')
%             disp([session_stopPhaseDir ' not found. Skipping ' sessionList{i_startSession} '...'])
%             continue
%         end
%         z_DiffmatName = ['region_z_stopSuccvFail_' sessionList{i_startSession} '.mat'];
%         z_DiffmatName = fullfile(session_stopPhaseDir, z_DiffmatName);
%         if ~exist(z_DiffmatName,'file'); continue; end
%         load(z_DiffmatName);
%         
%         f = z_STOPpower_metadata.f; t = z_STOPpower_metadata.t;
%         numEventTypes = length(z_STOPpower_metadata.eventList);
%         numFreqs  = length(f);
%         numSamps  = length(t);
%         
%         break;
%     end
%         

    session_stopPowerDir = fullfile(subject_stopPowerDir, sessionList{1});
    z_power_matName = ['power_stopSuccvFail_z_' STOPchannels{1}.name '_gabor.mat'];
    z_power_matName = fullfile(session_stopPowerDir, z_power_matName);
    load(z_power_matName);
    eventList = z_STOPpower_metadata.eventList;
    numEventTypes = length(eventList);
    f = z_STOPpower_metadata.f;t = z_STOPpower_metadata.t;
    f_to_plot = f;   % can change this later to be within a certain limit
    frange = [0,510];
    desired_plot_f_ticks = desired_freq_ticks(desired_freq_ticks > frange(1) & ...
                                          desired_freq_ticks < frange(2));
    freqTick_idx = zeros(1,length(desired_plot_f_ticks));
    for i_freqTick = 1 : length(desired_plot_f_ticks)
        freqTick_idx(i_freqTick) = find(abs(f_to_plot - desired_plot_f_ticks(i_freqTick)) == ...  
                                            min(abs(f_to_plot - desired_plot_f_ticks(i_freqTick))));
    end
    
    numFreqs = length(f);numSamps = length(t);
    
    region_z = zeros(num_allRegions, numEventTypes, numSamps, numFreqs);
    region_powerDiff = zeros(num_allRegions, numEventTypes, numSamps, numFreqs);
    numRegionSessions = zeros(1, num_allRegions);
    
    for iSession = 1 : numSessions
        disp(sprintf('Session %d out of %d', iSession, numSessions));
        
        session_stopPowerDir = fullfile(subject_stopPowerDir, sessionList{iSession});
        if ~exist(session_stopPowerDir, 'dir')
            disp([session_stopPowerDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end


        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        chList = extractChannels(cp, STOPchannels);
        sessionChannels = STOPchannels(chList);
        
        sessionRegions = getRegionsfromChannelDB(sessionChannels);
        numSessionRegions = length(sessionRegions);
        
        STOPpower_z_acrossSessions_metadata(iSession).sessionRegions = sessionRegions;
        STOPpower_z_acrossSessions_metadata(iSession).allRegions = allRegionList;
        STOPpower_z_acrossSessions_metadata(iSession).t = t;
        STOPpower_z_acrossSessions_metadata(iSession).f = f;
        STOPpower_z_acrossSessions_metadata(iSession).Fs = z_STOPpower_metadata.Fs;
        STOPpower_z_acrossSessions_metadata(iSession).eventList = z_STOPpower_metadata.eventList;
        STOPpower_z_acrossSessions_metadata(iSession).chName = z_STOPpower_metadata.chName;
        STOPpower_z_acrossSessions_metadata(iSession).twin = z_STOPpower_metadata.twin;
        
        numValidCh = zeros(1, numSessionRegions);
        for iRegion = 1 : numSessionRegions
            
            cp = initChanParams();
            cp.locationName = sessionRegions{iRegion};
            chList = extractChannels(cp, sessionChannels);
            regionChannels = sessionChannels(chList);
            numCh = length(regionChannels);
            allRegionIdx = find(strcmpi(sessionRegions{iRegion}, allRegionList));
        
            sessionRegion_z = NaN(numCh, numEventTypes, numSamps, numFreqs);
            sessionRegion_powerDiff = NaN(numCh, numEventTypes, numSamps, numFreqs);
            
            for iCh = 1 : numCh
            
                ch = sessionChannels{iCh};

                z_power_matName = ['power_stopSuccvFail_z_' ch.name '_gabor.mat'];
                z_power_matName = fullfile(session_stopPowerDir, z_power_matName);
                if ~exist(z_power_matName,'file'); continue; end
                
                numValidCh(iRegion) = numValidCh(iRegion) + 1;
                load(z_power_matName);
                
                sessionRegion_z(iCh, :, :, :) = zPowerDiff;
                sessionRegion_powerDiff(iCh, :, :, :) = realPowerDiff;
            end
            
            if numValidCh(iRegion) > 0
                region_z(allRegionIdx,:,:,:) = squeeze(region_z(allRegionIdx,:,:,:)) + squeeze(nanmean(sessionRegion_z,1));
                region_powerDiff(allRegionIdx,:,:,:) = squeeze(region_powerDiff(allRegionIdx,:,:,:)) + squeeze(nanmean(sessionRegion_powerDiff,1));
                numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
            end
        end
        
    end
    
    for i_allRegion = 1 : num_allRegions
        
        region_z(i_allRegion,:,:,:) = squeeze(region_z(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
        region_powerDiff(i_allRegion,:,:,:) = squeeze(region_powerDiff(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
        
%         numSessionRegions = length(STOPpower_z_acrossSessions_metadata(iSession).sessionRegions);
% 
%         for iSessionRegion = 1 : numSessionRegions
%             
%             allRegionIdx = find(strcmpi(STOPpower_z_acrossSessions_metadata.sessionRegions{iSessionRegion}, ...
%                                         allRegionList));
% %             mean_z = re_mean_z + 1i*im_mean_z;                        
%             region_z(allRegionIdx,:,:,:) = region_z(allRegionIdx,:,:,:) + ...
%                                            mean_z(iSessionRegion,:,:,:);
%             numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
%             
%         end

    end
    
    figProps.n = numEventTypes;
    figProps.colSpacing = ones(1, figProps.n) * 0.5;
    figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                 2 * sideMargins - ...
                                                 sum(figProps.colSpacing)) / figProps.n;

    zphase_acrossSessions_fig_saveName = [implantID '_zphase_stopSuccvFail_gabor'];
    zphase_acrossSessions_fig_saveName = fullfile(subject_stopPowerDir, zphase_acrossSessions_fig_saveName);
    zphase_acrossSessions_fig_saveName_zoom = [implantID '_zphase_stopSuccvFail_gabor_zoom'];
    zphase_acrossSessions_fig_saveName_zoom = fullfile(subject_stopPowerDir, zphase_acrossSessions_fig_saveName_zoom);
    
    power_z_acrossSessions_mat_saveName = [implantID '_power_z_stopSuccvFail_gabor.mat'];
    power_z_acrossSessions_mat_saveName = fullfile(subject_stopPowerDir, power_z_acrossSessions_mat_saveName);
        
    numRegionPlots = 0;
    numPages = 0;
    for iRegion = 1 : num_allRegions
%         cp = initChanParams();
%         cp.locationName = allRegionList{iRegion};            
%         chList = extractChannels(cp, sessionChannels);
%         if isempty(chList); continue; end
%         regionChannels = sessionChannels(chList);
%         numRegionChannels = length(regionChannels);

%         region_z(iRegion,:,:,:) = region_z(iRegion,:,:,:) / numRegionSessions(iRegion);

%         mrvDiff = squeeze(region_mrv(iRegion,1,:,:,:) - region_mrv(iRegion,2,:,:,:));
                
        numRegionPlots = numRegionPlots + 1;
        rowNum = rem(numRegionPlots, regions_per_page);
        if rowNum == 1
            [h_power_fig, h_power_axes] = createFigPanels5(figProps);
            [h_power_fig_zoom, h_power_axes_zoom] = createFigPanels5(figProps);
            page_regionList = allRegionList{iRegion};   %ch.name;
            page_numSessionList = num2str(numRegionSessions(iRegion));
            numPages = numPages + 1;
        else
            page_regionList = [page_regionList ', ' allRegionList{iRegion}];
            page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
        end
        if rowNum == 0; rowNum = regions_per_page; end

        for iEventType = 1 : numEventTypes
            axes(h_power_axes(rowNum, iEventType));

            f_to_plot = f;
            z_toPlot = squeeze(region_z(iRegion,iEventType, :, :))';
%             imagesc(t,1:length(f),z_toPlot);
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,z_toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', z_colorLim);
            set(gca,'yscale','log');
            
            if rowNum == 1
                title(z_STOPpower_metadata.eventList{iEventType});
            end
            if rowNum < regions_per_page
                set(gca,'xticklabel',[]);
            end
            if iEventType > 1
                set(gca,'yticklabel',[]);
            else
                ylabel('frequency (Hz)')
                set(gca,'ytick',round(f_to_plot(freqTick_idx)));
%                 set(gca,'yticklabel',round(f_to_plot(freqTick_idx)));
            end
            
            axes(h_power_axes_zoom(rowNum, iEventType));
            
            zoom_f_idx = (f>y_zoom_lim(1)) & (f<y_zoom_lim(2));
            z_toPlot = squeeze(region_z(iRegion,iEventType, :, zoom_f_idx))';
            f_to_plot = f(f>y_zoom_lim(1) & f<y_zoom_lim(2));
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,z_toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', z_colorLim,'ylim',y_zoom_lim);
            set(gca,'yscale','log');
            
            if rowNum == 1
                title(z_STOPpower_metadata.eventList{iEventType});
            end
            if rowNum < regions_per_page
                set(gca,'xticklabel',[]);
            end
            if iEventType > 1
                set(gca,'yticklabel',[]);
            else
                ylabel('frequency (Hz)')
                fticks = f(freqTick_idx);
                fticks = fticks((fticks > y_zoom_lim(1)) & (fticks < y_zoom_lim(2)));
%                 set(gca,'ytick',round(f_to_plot(freqTick_idx)));
                set(gca,'ytick',round(fticks));
            end
        end

        if (rowNum == regions_per_page || iRegion == num_allRegions)

            h_figAxes = createFigAxes(h_power_fig);
            h_figAxes_zoom = createFigAxes(h_power_fig_zoom);

            textStr{1} = ['difference in power z-scores (correct - failed), averaged across sessions'];
            textStr{2} = page_regionList;
            textStr{3} = ['number of sessions per region: ' page_numSessionList];
            textStr{4} = sprintf('color limits: %d to %d',z_colorLim(1),z_colorLim(2));
            
            axes(h_figAxes);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            svnm = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName,numPages);
            svnm_fig = sprintf('%s_%02d.fig',zphase_acrossSessions_fig_saveName,numPages);
            print(svnm,'-dpdf');
            
            axes(h_figAxes_zoom);
            svnm_zoom = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName_zoom,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm_zoom,'-dpdf');
%             if numPages == 1
                
%             else
%                 export_fig(zphase_acrossSessions_fig_saveName,'-pdf','-q101','-painters','-append');
%             end
            close(h_power_fig);
            close(h_power_fig_zoom);

        end

    end    % for iRegion...
    
    save(power_z_acrossSessions_mat_saveName, 'region_z', 'allRegionList', 'numRegionSessions', 'STOPpower_z_acrossSessions_metadata');
        
end    % for i_chDB...


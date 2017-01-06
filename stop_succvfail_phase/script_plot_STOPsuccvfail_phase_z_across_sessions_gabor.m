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
abs_colorLim = [0 2];
phaseDiff_colorLim = [0 1];
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
    
    for i_startSession = 1 : numSessions
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{i_startSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{i_startSession} '...'])
            continue
        end
        z_phase_matName = ['vecDiff_z_stopSuccvFail_' sessionList{i_startSession} '_gabors.mat'];
        z_phase_matName = fullfile(session_stopPhaseDir, z_phase_matName);
        if ~exist(z_phase_matName,'file'); continue; end
        load(z_phase_matName);
        
        f = STOPmrv_z_metadata.f; t = STOPmrv_z_metadata.t;
        eventList = STOPmrv_z_metadata.eventList;
        numEventTypes = length(eventList);
%         numFreqs  = length(f);
%         numSamps  = length(t);
        
        break;
    end
        

%     session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{1});
%     z_phase_matName = ['vecDiff_z_stopSuccvFail_' sessionList{1} '_gabors.mat'];
%     z_phase_matName = fullfile(session_stopPhaseDir, z_phase_matName);
%     load(z_phase_matName);
%     eventList = STOPmrv_z_metadata.eventList;
%     numEventTypes = length(eventList);
    f = STOPmrv_z_metadata.f;t = STOPmrv_z_metadata.t;
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
    region_mrv_diff = zeros(num_allRegions, numEventTypes, numSamps, numFreqs);
    region_mrv_abs_diff = zeros(num_allRegions, numEventTypes, numSamps, numFreqs);
    numRegionSessions = zeros(1, num_allRegions);
    
    for iSession = 1 : numSessions
        disp(sprintf('Session %d out of %d', iSession, numSessions));
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end


        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        chList = extractChannels(cp, STOPchannels);
        sessionChannels = STOPchannels(chList);
        
        sessionRegions = getRegionsfromChannelDB(sessionChannels);
        numSessionRegions = length(sessionRegions);
        
        z_phase_matName = ['vecDiff_z_stopSuccvFail_' sessionList{iSession} '_gabors.mat'];
        z_phase_matName = fullfile(session_stopPhaseDir, z_phase_matName);
        if ~exist(z_phase_matName,'file'); continue; end

        phase_matName = ['vecDiff_stopSuccvFail_' sessionList{iSession} '_gabors.mat'];
        phase_matName = fullfile(session_stopPhaseDir, phase_matName);
        if ~exist(phase_matName,'file'); continue; end
        load(phase_matName);
        load(z_phase_matName);
                
        STOPphase_z_acrossSessions_metadata(iSession).sessionRegions = sessionRegions;
        STOPphase_z_acrossSessions_metadata(iSession).allRegions = allRegionList;
        STOPphase_z_acrossSessions_metadata(iSession).t = t;
        STOPphase_z_acrossSessions_metadata(iSession).f = f;
        STOPphase_z_acrossSessions_metadata(iSession).Fs = STOPmrv_z_metadata.Fs;
        STOPphase_z_acrossSessions_metadata(iSession).eventList = STOPmrv_z_metadata.eventList;
        STOPphase_z_acrossSessions_metadata(iSession).chName = STOPmrv_z_metadata.chNames;
        STOPphase_z_acrossSessions_metadata(iSession).twin = STOPmrv_z_metadata.twin;
        
%         numValidCh = zeros(1, numSessionRegions);
        numValidCh = STOPmrv_z_metadata.num_ch_per_region;
        for iRegion = 1 : numSessionRegions
            
            cp = initChanParams();
            cp.locationName = sessionRegions{iRegion};
            chList = extractChannels(cp, sessionChannels);
            regionChannels = sessionChannels(chList);
            numCh = length(regionChannels);
            allRegionIdx = find(strcmpi(sessionRegions{iRegion}, allRegionList));
        
%             sessionRegion_z = NaN(numCh, numEventTypes, numSamps, numFreqs);
%             sessionRegion_powerDiff = NaN(numCh, numEventTypes, numSamps, numFreqs);
            
%             for iCh = 1 : numCh
%             
%                 ch = sessionChannels{iCh};
% 
%                 z_phase_matName = ['power_stopSuccvFail_z_' ch.name '_gabor.mat'];
%                 z_phase_matName = fullfile(session_stopPowerDir, z_phase_matName);
%                 if ~exist(z_phase_matName,'file'); continue; end
%                 
%                 numValidCh(iRegion) = numValidCh(iRegion) + 1;
%                 load(z_phase_matName);
%                 
%                 sessionRegion_z(iCh, :, :, :) = zPowerDiff;
%                 sessionRegion_powerDiff(iCh, :, :, :) = realPowerDiff;
%             end
            
            if numValidCh(iRegion) > 0
                region_z(allRegionIdx,:,:,:) = squeeze(region_z(allRegionIdx,:,:,:)) + squeeze(mean_mrv_z(iRegion,:,:,:));
                % CORRECT STOP - FAILED STOP
                mrv_diff = squeeze(re_mean_mrv(iRegion,1,:,:,:) - re_mean_mrv(iRegion,2,:,:,:) + ...
                                   1i * (im_mean_mrv(iRegion,1,:,:,:) - im_mean_mrv(iRegion,2,:,:,:)));
                region_mrv_diff(allRegionIdx,:,:,:) = squeeze(region_mrv_diff(allRegionIdx,:,:,:)) + mrv_diff;
%                 mrvDiffAmp = abs(squeeze(re_mean_mrv(iRegion,:,:,:)) + 1i * squeeze(im_mean_mrv(iRegion,:,:,:)));
                region_mrv_abs_diff(allRegionIdx,:,:,:) = squeeze(region_mrv_abs_diff(allRegionIdx,:,:,:)) + abs(mrv_diff);
                numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
            end
        end
        
    end
    
    for i_allRegion = 1 : num_allRegions
        
        region_z(i_allRegion,:,:,:) = squeeze(region_z(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
        region_mrv_diff(i_allRegion,:,:,:) = squeeze(region_mrv_diff(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
        region_mrv_abs_diff(i_allRegion,:,:,:) = squeeze(region_mrv_abs_diff(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
        
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

    zphase_acrossSessions_fig_saveName = [implantID '_zphase_stopSuccvFail_gabor.pdf'];
    zphase_acrossSessions_fig_saveName = fullfile(subject_stopPhaseDir, zphase_acrossSessions_fig_saveName);
    zphase_acrossSessions_fig_saveName_zoom = [implantID '_zphase_stopSuccvFail_gabor_zoom'];
    zphase_acrossSessions_fig_saveName_zoom = fullfile(subject_stopPhaseDir, zphase_acrossSessions_fig_saveName_zoom);
    abs_phaseDiff_acrossSessions_fig_saveName = [implantID '_abs_mrvDiff_stopSuccvFail_gabor'];
    abs_phaseDiff_acrossSessions_fig_saveName = fullfile(subject_stopPhaseDir, abs_phaseDiff_acrossSessions_fig_saveName);
    abs_phaseDiff_acrossSessions_fig_saveName_zoom = [implantID '_abs_mrvDiff_stopSuccvFail_gabor_zoom'];
    abs_phaseDiff_acrossSessions_fig_saveName_zoom = fullfile(subject_stopPhaseDir, abs_phaseDiff_acrossSessions_fig_saveName_zoom);
    phaseDiff_acrossSessions_fig_saveName = [implantID '_meanMRVdiff_stopSuccvFail_gabor'];
    phaseDiff_acrossSessions_fig_saveName = fullfile(subject_stopPhaseDir, phaseDiff_acrossSessions_fig_saveName);
    phaseDiff_acrossSessions_fig_saveName_zoom = [implantID '_meanMRVdiff_stopSuccvFail_gabor_zoom'];
    phaseDiff_acrossSessions_fig_saveName_zoom = fullfile(subject_stopPhaseDir, phaseDiff_acrossSessions_fig_saveName_zoom);
    
    fig_saveName = [implantID '_stopSuccvFail_gabor_phases'];
    fig_saveName = fullfile(subject_stopPhaseDir, fig_saveName);
    
    phase_z_acrossSessions_mat_saveName = [implantID '_phaseDiff_z_stopSuccvFail_gabor.mat'];
    phase_z_acrossSessions_mat_saveName = fullfile(subject_stopPhaseDir, phase_z_acrossSessions_mat_saveName);
        
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
            [h_phase_z_fig, h_phase_z_axes] = createFigPanels5(figProps);
            [h_phase_z_fig_zoom, h_phase_z_axes_zoom] = createFigPanels5(figProps);
            
            [h_phase_amp_fig, h_phase_amp_axes] = createFigPanels5(figProps);
            [h_phase_amp_fig_zoom, h_phase_amp_axes_zoom] = createFigPanels5(figProps);
            
            [h_phaseDiff_fig, h_phaseDiff_axes] = createFigPanels5(figProps);
            [h_phaseDiff_fig_zoom, h_phaseDiff_axes_zoom] = createFigPanels5(figProps);
            
            page_regionList = allRegionList{iRegion};   %ch.name;
            page_numSessionList = num2str(numRegionSessions(iRegion));
            numPages = numPages + 1;
        else
            page_regionList = [page_regionList ', ' allRegionList{iRegion}];
            page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
        end
        if rowNum == 0; rowNum = regions_per_page; end

        for iEventType = 1 : numEventTypes
            axes(h_phase_z_axes(rowNum, iEventType));

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
            
            axes(h_phase_z_axes_zoom(rowNum, iEventType));
            
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
            
            
            
            % DIFFERENCES IN MRV ABSOLUTE VALUE
            axes(h_phase_amp_axes(rowNum, iEventType));

            f_to_plot = f;
            toPlot = squeeze(region_mrv_abs_diff(iRegion,iEventType, :, :))';
%             imagesc(t,1:length(f),z_toPlot);
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', abs_colorLim);
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
        
            % DIFFERENCES IN MRV ABSOLUTE VALUE, ZOOMED
            axes(h_phase_amp_axes_zoom(rowNum, iEventType));

            f_to_plot = f(f>y_zoom_lim(1) & f<y_zoom_lim(2));
            toPlot = squeeze(region_mrv_abs_diff(iRegion,iEventType, :, zoom_f_idx))';
%             imagesc(t,1:length(f),z_toPlot);
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', abs_colorLim);
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
            
            
            % DIFFERENCES IN PHASE (SUBTRACTING COMPLEX VECTORS AND 
            % AVERAGING INSTEAD OF THE MAGNITUDE OF EACH MRV)
            axes(h_phaseDiff_axes(rowNum, iEventType));

            f_to_plot = f;
            toPlot = abs(squeeze(region_mrv_diff(iRegion,iEventType, :, :)))';
%             imagesc(t,1:length(f),z_toPlot);
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', phaseDiff_colorLim);
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
            

            axes(h_phaseDiff_axes_zoom(rowNum, iEventType));

            f_to_plot = f(f>y_zoom_lim(1) & f<y_zoom_lim(2));
            toPlot = abs(squeeze(region_mrv_diff(iRegion,iEventType, :, zoom_f_idx)))';
%             imagesc(t,1:length(f),z_toPlot);
%             imagesc(t,f_to_plot,z_toPlot);
            h_pcolor = pcolor(t,f_to_plot,toPlot);
            h_pcolor.EdgeColor = 'none';
            set(gca,'ydir','normal','clim', abs_colorLim);
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

            h_figAxes = createFigAxes(h_phase_z_fig);
            h_figAxes_zoom = createFigAxes(h_phase_z_fig_zoom);
            h_figAxes_abs = createFigAxes(h_phase_amp_fig);
            h_figAxes_abs_zoom = createFigAxes(h_phase_amp_fig_zoom);
            h_figAxes_phaseDiff = createFigAxes(h_phaseDiff_fig);
            h_figAxes_phaseDiff_zoom = createFigAxes(h_phaseDiff_fig_zoom);

            textStr{1} = implantID;
            textStr{2} = ['difference in z-score phase (correct - failed), averaged across sessions'];
            textStr{3} = page_regionList;
            textStr{4} = ['number of sessions per region: ' page_numSessionList];
            textStr{5} = sprintf('color limits: %d to %d',z_colorLim(1),z_colorLim(2));
            
            axes(h_figAxes);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            svnm = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName,numPages);
            print(svnm,'-dpdf');
            
            axes(h_figAxes_zoom);
            svnm_zoom = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName_zoom,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm_zoom,'-dpdf');
            
            
            textStr{2} = ['absolute difference mrv values (correct - failed), averaged across sessions'];
            textStr{5} = sprintf('color limits: %d to %d',abs_colorLim(1),abs_colorLim(2));
            axes(h_figAxes_abs);
            svnm = sprintf('%s_%02d.pdf',abs_phaseDiff_acrossSessions_fig_saveName,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm,'-dpdf');
            
            axes(h_figAxes_abs_zoom);
            svnm_zoom = sprintf('%s_%02d.pdf',abs_phaseDiff_acrossSessions_fig_saveName_zoom,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm_zoom,'-dpdf');
            
            
            textStr{2} = ['difference between mrv vectors (correct - failed), averaged across sessions'];
            textStr{5} = sprintf('color limits: %d to %d',phaseDiff_colorLim(1),phaseDiff_colorLim(2));
            axes(h_figAxes_phaseDiff);
            svnm = sprintf('%s_%02d.pdf',phaseDiff_acrossSessions_fig_saveName,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm,'-dpdf');
            
            axes(h_figAxes_phaseDiff_zoom);
            svnm_zoom = sprintf('%s_%02d.pdf',phaseDiff_acrossSessions_fig_saveName_zoom,numPages);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            print(svnm_zoom,'-dpdf');
            
            figName = sprintf('%s_%02d.fig',fig_saveName,numPages);
            savefig([h_phase_z_fig, h_phase_amp_fig, h_phaseDiff_fig], figName, 'compact');
%             if numPages == 1
                
%             else
%                 export_fig(zphase_acrossSessions_fig_saveName,'-pdf','-q101','-painters','-append');
%             end
            close(h_phase_z_fig);
            close(h_phase_z_fig_zoom);
            close(h_phase_amp_fig);
            close(h_phase_amp_fig_zoom);
            close(h_phaseDiff_fig);
            close(h_phaseDiff_fig_zoom);

        end

    end    % for iRegion...
    
    save(phase_z_acrossSessions_mat_saveName, 'region_z', 'region_mrv_diff', 'region_mrv_abs_diff', 'allRegionList', 'numRegionSessions', 'STOPpower_z_acrossSessions_metadata');
        
end    % for i_chDB...


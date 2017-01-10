% script_plotRT_phaseCorrelations_gabor_across_sessions

% script to go through all STOP sessions and see if there are consistent
% phase differences between STOP-success and STOP-failure trials

% chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
% phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';

chDB_directory         = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
phaseRTcorr_directory  = '/Volumes/Tbolt_02/stop-signal reanalysis/phase_RT_correlations_gabors';
RTcorr_plots_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/phase_RT correlation gabor plots';

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

y_zoom_lim = [1 100];   % for more focused plot
y_full_plot_lim = [0.1 502];

regions_per_page = 5;
colorLim = [0 1];

% desired_freq_ticks = [4,8,16,32,64,128,256,500];
desired_freq_ticks = [4,16,64,128,500];

for i_chDB = 1 : 4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_gabor_plots']);
    if ~exist(subject_RTcorr_plots_directory, 'dir')
        disp([subject_RTcorr_plots_directory ' not found. Skipping ' implantID '...'])
        continue
    end

    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
        
    % COMMENT OUT THESE LINES TO LOOK AT ALL REGIONS
    cp = initChanParams();
    cp.locationSubClass = ROI_list;    
    chList = extractChannels( cp, channels );
    channels = channels(chList);
            
            
            
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
    channels = excludeChannels(cp, channels);
    
    cp = initChanParams();
    cp.locationSubClass = {'e2', 'e3', 'e02','e03'};
    channels = excludeChannels(cp, channels);
    if isempty(channels);continue;end

    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    allRegionList = getSubclassesfromChannelDB( channels );
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
%         f = z_STOPphase_metadata.f; t = z_STOPphase_metadata.t;
%         numEventTypes = length(z_STOPphase_metadata.eventList);
%         numFreqs  = length(f);
%         numSamps  = length(t);
%         
%         break;
%     end
%         

    session_RTcorr_plots_directory = fullfile(subject_RTcorr_plots_directory, [sessionList{1}, '_RTcorr_gabor_plots']);
    session_RTphase_matName = [sessionList{1} '_phase_RTcorr_gabor'];
    session_RTphase_matName = fullfile(session_RTcorr_plots_directory, session_RTphase_matName);
    load(session_RTphase_matName);
    eventList = session_phaseRTcorr_metadata.eventList;
    numEventTypes = length(eventList);
    f = session_phaseRTcorr_metadata.freqs;t = session_phaseRTcorr_metadata.t;
    f_to_plot = f;   % can change this later to be within a certain limit
    frange = [0,510];
    desired_plot_f_ticks = desired_freq_ticks(desired_freq_ticks > frange(1) & ...
                                          desired_freq_ticks < frange(2));
    freqTick_idx = zeros(1,length(desired_plot_f_ticks));
    for i_freqTick = 1 : length(desired_plot_f_ticks)
        freqTick_idx(i_freqTick) = find(abs(f_to_plot - desired_plot_f_ticks(i_freqTick)) == ...  
                                            min(abs(f_to_plot - desired_plot_f_ticks(i_freqTick))));
    end
    
    % set up plots
    figProps.width  = 11 * 2.54;
    figProps.height = 8.5 * 2.54;

    figProps.m = regions_per_page;    % number of rows
    figProps.n = numEventTypes;

    sideMargins = 2; botMargin = 2.54;

    figProps.rowSpacing = ones(1, figProps.m) * 0.5;
    figProps.colSpacing = ones(1, figProps.n) * 0.5;

    figProps.topMargin  = 5;

    figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                 2 * sideMargins - ...
                                                 sum(figProps.colSpacing)) / figProps.n;
    figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                  figProps.topMargin - botMargin - ...
                                                  sum(figProps.rowSpacing)) / figProps.m;
                                          
                                          
    numFreqs = length(f);numSamps = length(t);
    
    region_RTphase_mean = zeros(num_allRegions, numEventTypes, numFreqs, numSamps);
    numRegionSessions = zeros(1, num_allRegions);
    
    for iSession = 1 : numSessions
        fprintf('Session %d out of %d\n', iSession, numSessions);
        
        session_RTcorr_plots_directory = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession}, '_RTcorr_gabor_plots']);
        if ~exist(session_RTcorr_plots_directory, 'dir')
            disp([session_RTcorr_plots_directory ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end
        
        session_RTphase_matName = [sessionList{iSession} '_phase_RTcorr_gabor.mat'];
        session_RTphase_matName = fullfile(session_RTcorr_plots_directory, session_RTphase_matName);
        
        if ~exist(session_RTphase_matName,'file'); continue; end
        load(session_RTphase_matName);

%         cp = initChanParams();
%         cp.session = sessionList{iSession};
%         chList = extractChannels(cp, channels);
%         sessionChannels = channels(chList);
        
        sessionRegions = session_phaseRTcorr_metadata.regionList;
        numSessionRegions = length(sessionRegions);
        
        phaseRTcorr_acrossSessions_metadata(iSession).sessionRegions = sessionRegions;
        phaseRTcorr_acrossSessions_metadata(iSession).allRegions = allRegionList;
        phaseRTcorr_acrossSessions_metadata(iSession).sessionRegions = sessionRegions;
        phaseRTcorr_acrossSessions_metadata(iSession).t = t;
        phaseRTcorr_acrossSessions_metadata(iSession).f = f;
        phaseRTcorr_acrossSessions_metadata(iSession).Fs = session_phaseRTcorr_metadata.Fs;
        phaseRTcorr_acrossSessions_metadata(iSession).eventList = session_phaseRTcorr_metadata.eventList;
        phaseRTcorr_acrossSessions_metadata(iSession).twin = session_phaseRTcorr_metadata.twin;
        phaseRTcorr_acrossSessions_metadata(iSession).numChPerRegion = session_phaseRTcorr_metadata.numChPerRegion;
        
        for iRegion = 1 : numSessionRegions
            
            allRegionIdx = find(strcmpi(sessionRegions{iRegion}, allRegionList));
            
            if isempty(allRegionIdx); continue; end
        
            region_RTphase_mean(allRegionIdx, :, :, :) = squeeze(region_RTphase_mean(allRegionIdx, :, :, :)) + ...
                squeeze(mean_phaseRTcorr(iRegion,:,:,:));
            
            numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
            
        end
        
    end    % iSession
    
    for iRegion = 1 : num_allRegions
        if numRegionSessions(iRegion) > 0
            region_RTphase_mean(iRegion,:,:,:) = region_RTphase_mean(iRegion,:,:,:) / numRegionSessions(iRegion);
        end
    end
    
    for iSession = 1 : numSessions
        phaseRTcorr_acrossSessions_metadata(iSession).numSessionsPerRegion = numRegionSessions;
    end
    
    phaseRTcorr_summary_name = [implantID '_phase_RTcorr_gabor'];
    phaseRTcorr_summary_name = fullfile(subject_RTcorr_plots_directory,phaseRTcorr_summary_name);
    
    save(phaseRTcorr_summary_name, 'region_RTphase_mean', 'phaseRTcorr_acrossSessions_metadata');
    
    % now make the summary plots
    phaseRT_summary_plot_name = [implantID '_phase_RTcorr_gabor_summary'];
    
    numPages = 0;
    for i_allRegion = 1 : num_allRegions
        rowNum = rem(i_allRegion, figProps.m);
        if rowNum == 0;rowNum = figProps.m;end
        
        if rowNum == 1
            [h_phaseCorr_fig, h_phaseCorr_axes] = createFigPanels5(figProps);
            page_regionList = allRegionList{i_allRegion};
            page_sessionsPerRegion = num2str(numRegionSessions(i_allRegion));
            numPages = numPages + 1;
        else
            page_regionList = [page_regionList ', ' allRegionList{i_allRegion}];
            page_sessionsPerRegion = [page_sessionsPerRegion ', ' num2str(numRegionSessions(i_allRegion))];
        end
        
        for iEvent = 1 : numEventTypes
            axes(h_phaseCorr_axes(rowNum, iEvent));

            toPlot = squeeze(region_RTphase_mean(i_allRegion,iEvent,:,:));
            h_pcolor = pcolor(t, f, toPlot);
            h_pcolor.EdgeColor = 'none';                                            % WORKING HERE...

            set(gca,'ydir','normal');
            set(gca,'yscale','log');
            set(gca,'clim',colorLim);
        
            if rowNum == 1
                title([eventList{iEvent}]);
            end

            if rowNum < figProps.m
                set(gca,'xticklabel',[]);
            end
            if iEvent > 1
                set(gca,'yticklabel',[]);
            else
                ylabel('Frequency (Hz)');
                set(gca,'ytick',round(f_to_plot(freqTick_idx)));
            end

        end
        
        if rem(i_allRegion, figProps.m) == 0 || i_allRegion == num_allRegions
            h_figAxes = createFigAxes(h_phaseCorr_fig);

            textStr{1} = ['Circular correlation coefficient between phase and RT averaged across sessions'];
            textStr{2} = page_regionList;
            textStr{3} = sprintf('valid sessions for each region: %s',page_sessionsPerRegion);
%                 textStr{4} = ['number of channels per region: ' page_numChList];
            textStr{4} = sprintf('color limits: %d to %d',colorLim(1),colorLim(2));

            cur_PDFname = sprintf('%s_%02d.pdf',phaseRT_summary_plot_name,numPages);
            cur_PDFname = fullfile(subject_RTcorr_plots_directory, cur_PDFname);
            cur_figName = sprintf('%s_%02d.fig',phaseRT_summary_plot_name,numPages);
            cur_figName = fullfile(subject_RTcorr_plots_directory, cur_figName);

            axes(h_figAxes);
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

            print(cur_PDFname, '-dpdf');
            savefig(h_phaseCorr_fig,cur_figName,'compact');

            close(h_phaseCorr_fig);
        end
        
    end

end
            
%             if numValidCh(iRegion) > 0
%                 region_z(allRegionIdx,:,:,:) = squeeze(region_z(allRegionIdx,:,:,:)) + squeeze(nanmean(sessionRegion_z,1));
%                 region_phaseDiff(allRegionIdx,:,:,:) = squeeze(region_phaseDiff(allRegionIdx,:,:,:)) + squeeze(nanmean(sessionRegion_phaseDiff,1));
%                 numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
%             end
%         end
%         
%     end
%     
%     for i_allRegion = 1 : num_allRegions
%         
%         region_z(i_allRegion,:,:,:) = squeeze(region_z(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
%         region_phaseDiff(i_allRegion,:,:,:) = squeeze(region_phaseDiff(i_allRegion,:,:,:)) / numRegionSessions(i_allRegion);
%         
% %         numSessionRegions = length(STOPphase_z_acrossSessions_metadata(iSession).sessionRegions);
% % 
% %         for iSessionRegion = 1 : numSessionRegions
% %             
% %             allRegionIdx = find(strcmpi(STOPphase_z_acrossSessions_metadata.sessionRegions{iSessionRegion}, ...
% %                                         allRegionList));
% % %             mean_z = re_mean_z + 1i*im_mean_z;                        
% %             region_z(allRegionIdx,:,:,:) = region_z(allRegionIdx,:,:,:) + ...
% %                                            mean_z(iSessionRegion,:,:,:);
% %             numRegionSessions(allRegionIdx) = numRegionSessions(allRegionIdx) + 1;
% %             
% %         end
% 
%     end
%     
%     figProps.n = numEventTypes;
%     figProps.colSpacing = ones(1, figProps.n) * 0.5;
%     figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
%                                                  2 * sideMargins - ...
%                                                  sum(figProps.colSpacing)) / figProps.n;
% 
%     zphase_acrossSessions_fig_saveName = [implantID '_zphase_stopSuccvFail_gabor'];
%     zphase_acrossSessions_fig_saveName = fullfile(subject_stopphaseDir, zphase_acrossSessions_fig_saveName);
%     zphase_acrossSessions_fig_saveName_zoom = [implantID '_zphase_stopSuccvFail_gabor_zoom'];
%     zphase_acrossSessions_fig_saveName_zoom = fullfile(subject_stopphaseDir, zphase_acrossSessions_fig_saveName_zoom);
%     
%     phase_z_acrossSessions_mat_saveName = [implantID '_phase_z_stopSuccvFail_gabor.mat'];
%     phase_z_acrossSessions_mat_saveName = fullfile(subject_stopphaseDir, phase_z_acrossSessions_mat_saveName);
%         
%     numRegionPlots = 0;
%     numPages = 0;
%     for iRegion = 1 : num_allRegions
% %         cp = initChanParams();
% %         cp.locationName = allRegionList{iRegion};            
% %         chList = extractChannels(cp, sessionChannels);
% %         if isempty(chList); continue; end
% %         regionChannels = sessionChannels(chList);
% %         numRegionChannels = length(regionChannels);
% 
% %         region_z(iRegion,:,:,:) = region_z(iRegion,:,:,:) / numRegionSessions(iRegion);
% 
% %         mrvDiff = squeeze(region_mrv(iRegion,1,:,:,:) - region_mrv(iRegion,2,:,:,:));
%                 
%         numRegionPlots = numRegionPlots + 1;
%         rowNum = rem(numRegionPlots, regions_per_page);
%         if rowNum == 1
%             [h_phase_fig, h_phase_axes] = createFigPanels5(figProps);
%             [h_phase_fig_zoom, h_phase_axes_zoom] = createFigPanels5(figProps);
%             page_regionList = allRegionList{iRegion};   %ch.name;
%             page_numSessionList = num2str(numRegionSessions(iRegion));
%             numPages = numPages + 1;
%         else
%             page_regionList = [page_regionList ', ' allRegionList{iRegion}];
%             page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
%         end
%         if rowNum == 0; rowNum = regions_per_page; end
% 
%         for iEventType = 1 : numEventTypes
%             axes(h_phase_axes(rowNum, iEventType));
% 
%             f_to_plot = f;
%             z_toPlot = squeeze(region_z(iRegion,iEventType, :, :))';
% %             imagesc(t,1:length(f),z_toPlot);
% %             imagesc(t,f_to_plot,z_toPlot);
%             h_pcolor = pcolor(t,f_to_plot,z_toPlot);
%             h_pcolor.EdgeColor = 'none';
%             set(gca,'ydir','normal','clim', z_colorLim);
%             set(gca,'yscale','log');
%             
%             if rowNum == 1
%                 title(z_STOPphase_metadata.eventList{iEventType});
%             end
%             if rowNum < regions_per_page
%                 set(gca,'xticklabel',[]);
%             end
%             if iEventType > 1
%                 set(gca,'yticklabel',[]);
%             else
%                 ylabel('frequency (Hz)')
%                 set(gca,'ytick',round(f_to_plot(freqTick_idx)));
% %                 set(gca,'yticklabel',round(f_to_plot(freqTick_idx)));
%             end
%             
%             axes(h_phase_axes_zoom(rowNum, iEventType));
%             
%             zoom_f_idx = (f>y_zoom_lim(1)) & (f<y_zoom_lim(2));
%             z_toPlot = squeeze(region_z(iRegion,iEventType, :, zoom_f_idx))';
%             f_to_plot = f(f>y_zoom_lim(1) & f<y_zoom_lim(2));
% %             imagesc(t,f_to_plot,z_toPlot);
%             h_pcolor = pcolor(t,f_to_plot,z_toPlot);
%             h_pcolor.EdgeColor = 'none';
%             set(gca,'ydir','normal','clim', z_colorLim,'ylim',y_zoom_lim);
%             set(gca,'yscale','log');
%             
%             if rowNum == 1
%                 title(z_STOPphase_metadata.eventList{iEventType});
%             end
%             if rowNum < regions_per_page
%                 set(gca,'xticklabel',[]);
%             end
%             if iEventType > 1
%                 set(gca,'yticklabel',[]);
%             else
%                 ylabel('frequency (Hz)')
%                 fticks = f(freqTick_idx);
%                 fticks = fticks((fticks > y_zoom_lim(1)) & (fticks < y_zoom_lim(2)));
% %                 set(gca,'ytick',round(f_to_plot(freqTick_idx)));
%                 set(gca,'ytick',round(fticks));
%             end
%         end
% 
%         if (rowNum == regions_per_page || iRegion == num_allRegions)
% 
%             h_figAxes = createFigAxes(h_phase_fig);
%             h_figAxes_zoom = createFigAxes(h_phase_fig_zoom);
% 
%             textStr{1} = ['difference in phase z-scores (correct - failed), averaged across sessions'];
%             textStr{2} = page_regionList;
%             textStr{3} = ['number of sessions per region: ' page_numSessionList];
%             textStr{4} = sprintf('color limits: %d to %d',z_colorLim(1),z_colorLim(2));
%             
%             axes(h_figAxes);
%             text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
%             svnm = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName,numPages);
%             svnm_fig = sprintf('%s_%02d.fig',zphase_acrossSessions_fig_saveName,numPages);
%             print(svnm,'-dpdf');
%             
%             axes(h_figAxes_zoom);
%             svnm_zoom = sprintf('%s_%02d.pdf',zphase_acrossSessions_fig_saveName_zoom,numPages);
%             text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
%             print(svnm_zoom,'-dpdf');
% %             if numPages == 1
%                 
% %             else
% %                 export_fig(zphase_acrossSessions_fig_saveName,'-pdf','-q101','-painters','-append');
% %             end
%             close(h_phase_fig);
%             close(h_phase_fig_zoom);
% 
%         end
% 
%     end    % for iRegion...
%     
%     save(phase_z_acrossSessions_mat_saveName, 'region_z', 'allRegionList', 'numRegionSessions', 'STOPphase_z_acrossSessions_metadata');
%         
% end    % for i_chDB...
% 

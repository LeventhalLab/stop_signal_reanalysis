% script_plot_phaseRTcorr_gabor_summaries

% chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% phaseRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabors';
% RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabor_plots';

chDB_directory    = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
phase_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

colorLim = [0, 1];

desired_freq_ticks = [4,8,16,32,64,128,256,500];

ROI_list = {'eegorb','cpu','gp','stn','snr'};
frange = [0,10];
pdfSuffix = sprintf('_phaseRTcorr_gabors_%d_%d_Hz',frange(1),frange(2));

for i_chDB = 1 : 4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
%     if ~exist(chDB_list{i_chDB}, 'var')
%         chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
%         disp(['loading ' chDB_file]);
%         load( chDB_file );
%     end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));

    subject_phaseRTcorr_directory = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr_gabors']);
    if ~exist(subject_phaseRTcorr_directory, 'dir')
        disp([subject_phaseRTcorr_directory ' not found.']);
        continue;
    end

    regionSummaryMatName = [implantID '_phaseRTcorr_gabors_across_sessions.mat'];
    regionSummaryMatName = fullfile(subject_phaseRTcorr_directory, regionSummaryMatName);
    
    if ~exist(regionSummaryMatName,'file');
        continue;
    end
    
    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_gabor_RTcorr_plots']);
    if ~exist(subject_RTcorr_plots_directory, 'dir')
        mkdir(subject_RTcorr_plots_directory);
    end
    
    load(regionSummaryMatName);
    eventList = region_phase_RTcorr_metadata.eventList;
    numEvents = length(eventList);
    twin = region_phase_RTcorr_metadata.twin;
    t = region_phase_RTcorr_metadata.t;
    numValidSessions = region_phase_RTcorr_metadata.sessions_per_region;
    
    regionList = region_phase_RTcorr_metadata.regionList;
    validRegions = numValidSessions > 0;
    f = region_phase_RTcorr_metadata.f;
    f_lim_idx = zeros(1,2);
    f_lim_idx(1) = find(f > frange(1),1,'first');
    f_lim_idx(2) = find(f < frange(2),1,'last');
    f_to_plot = f(f_lim_idx(1):f_lim_idx(2));
    xticks = [twin(1),0,twin(2)];
    
    desired_plot_f_ticks = desired_freq_ticks(desired_freq_ticks > frange(1) & ...
                                              desired_freq_ticks < frange(2));
    freqTick_idx = zeros(1,length(desired_plot_f_ticks));
    for i_freqTick = 1 : length(desired_plot_f_ticks)
        freqTick_idx(i_freqTick) = find(abs(f_to_plot - desired_plot_f_ticks(i_freqTick)) == ...  
                                            min(abs(f_to_plot - desired_plot_f_ticks(i_freqTick))));
    end
    
    region_phaseRTcorr = mean_phaseRTcorr(validRegions,:,:,:);
    validRegionNames = cell(1,1);
    numValidRegions = 0;
    for ii = 1 : length(validRegions)
        if validRegions(ii)
            numValidRegions = numValidRegions + 1;
            validRegionNames{numValidRegions} = regionList{ii};
        end
    end
    
    figProps.width  = 11 * 2.54;
    figProps.height = 8.5 * 2.54;

    figProps.m = channels_per_page; figProps.n = numEvents;   % number of rows and columns, respectively

    sideMargins = 1.5; botMargin = 2.54;

    figProps.colSpacing = ones(1, figProps.n) * 0.5;
    figProps.rowSpacing = ones(1, figProps.m) * 0.5;
    figProps.topMargin  = 5;

    figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                 2 * sideMargins - ...
                                                 sum(figProps.colSpacing)) / figProps.n;
    figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                  figProps.topMargin - botMargin - ...
                                                  sum(figProps.rowSpacing)) / figProps.m;
    pdfBase = sprintf('%s%s',implantID,pdfSuffix);

    testName = sprintf('%s_01.pdf',pdfBase);
    testName = fullfile(subject_RTcorr_plots_directory, testName);
    if exist(testName,'file'); continue; end
    numPagesForRegions = 0;
    
    for iRegion = 1 : length(ROI_list)
        
        regionIdx = strcmpi(ROI_list{iRegion},validRegionNames);
        full_regionList_idx = strcmpi(ROI_list{iRegion}, region_phase_RTcorr_metadata.regionList);
        
        rowNum = rem(iRegion, channels_per_page);
        if rowNum == 0; rowNum = channels_per_page; end
        
        if rowNum == 1
            numPagesForRegions = numPagesForRegions + 1;
            [h_fig, h_axes] = createFigPanels5(figProps);
            page_regionList = ROI_list{iRegion};
            page_chPerRegion = num2str(numValidSessions(iRegion));
            pageName = sprintf('%s_%02d.pdf',pdfBase,numPagesForRegions);
            PDFname = fullfile(subject_RTcorr_plots_directory, pageName);
        else
            page_regionList = [page_regionList ', ' ROI_list{iRegion}];
            if any(regionIdx)
                page_chPerRegion = [page_chPerRegion ', ' num2str(numValidSessions(full_regionList_idx))];
            else
                page_chPerRegion = [page_chPerRegion ', 0'];
            end
        end
        if any(regionIdx)
        
            for iEventType = 1 : numEvents

                axes(h_axes(iRegion, iEventType));
                toPlot = squeeze(region_phaseRTcorr(regionIdx, iEventType, f_lim_idx(1):f_lim_idx(2), :));

                imagesc(t, f_to_plot, toPlot);
                set(gca,'ydir','normal');
                set(gca,'clim',colorLim);
                if rowNum == 1
                    title(eventList{iEventType});
                end

                if rowNum < figProps.m
                    set(gca,'xticklabel',{});
                else
                    xlabel('time (s)')
                    set(gca,'xtick',xticks);
                end
                if iEventType > 1
                    set(gca,'yticklabel',{});
                else
                    ylabel('frequency (Hz)')
                    set(gca,'ytick',round(f_to_plot(freqTick_idx)));
                end

            end    % for iEventType
            
        end
        
        if rem(iRegion, figProps.m) == 0 || iRegion == numRegions
            h_figAxes = createFigAxes(h_fig);
            axes(h_figAxes);

            textStr{1} = [implantID ', phase-RT correlations (circstat), averages across sessions'];
            textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
            textStr{3} = page_regionList;
            textStr{4} = ['number of sessions per region: ' num2str(page_chPerRegion)];
            textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

%             if numPagesForRegions == 1
%                 export_fig(PDFname, '-pdf', '-q101', '-painters');
%             else
%                 export_fig(PDFname, '-pdf', '-q101', '-painters', '-append');
%             end
            print(PDFname, '-dpdf');
            close(h_fig);
        end
        
    end    % for iRegion...
    
    
end    % for i_chDB...
    
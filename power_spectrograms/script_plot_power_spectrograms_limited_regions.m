% script_plot_power_spectrograms

% UPDATED 01-12-2015

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

powerSpectrogramDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectrograms';
% phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

% eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};
% numEvents = length(eventList);


% eventList{1} = {'noseCenterIn'};
% eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
% eventList{3} = eventList{2};
% eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
% eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
% eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
% eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

% numROI = zeros(1, length(ROI_list));
% numValidRegions = length(ROI_list);

ROI_list = {'eegorb','cpu','gp','stn','snr'};                               % regions of interest

numTrialTypes = length(trialTypeList);

% eventtWin(1,:) = [-1 2];   % for analysis of all trials
% % analysisWin(1) = 3;
% % stepSize(1)    = 3;
% for iTrialType = 2 : numTrialTypes
%     eventtWin(iTrialType,:) = [-1 1];
% %     analysisWin(iTrialType) = 0.1;
% %     stepSize(iTrialType)    = 0.05;
% end
regions_per_page = 5;
log_clim = [-13 -2];
power_clim = [0, 1e-2];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;
figProps.m = regions_per_page;    % number of rows
sideMargins = 2; botMargin = 2.54;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
% NOW NEED TO AVERAGE ACROSS ANIMALS
                                          
for i_chDB = 1 : 4%length(chDB_list)
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    
    subject_powerSpectDir = fullfile(powerSpectrogramDir, [implantID '_powerSpectrograms']);
    if ~exist(subject_powerSpectDir, 'dir')
        continue;
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
%     channels = eval( chDB_info.name );
%     
%     cp = initChanParams();
%     cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
%     channels = excludeChannels(cp, channels);
    
    if i_chDB == 1
        startTrialType = 4;
    else
        startTrialType = 1;
    end
    for iTrialType = startTrialType : length(trialTypeList)
        if iTrialType == 2 || iTrialType == 3   % take this out once calculations updated for nose-side-out event
            continue;
        end
        
        trialType = trialTypeList{iTrialType}

%         powerSpect_metadata.eventList = eventList{iTrialType};
%         powerSpect_metadata.trialType = trialType;
%         powerSpect_metadata.eventtWin = twin;
        
        regionSummaryMatName = [implantID '_' trialType '_powerSpect_across_sessions.mat'];
        regionSummaryMatName = fullfile(subject_powerSpectDir, regionSummaryMatName);
        load(regionSummaryMatName);
        
        % mean_spect is a 4-D array: region x event x frequency x time
        eventList = region_power_spectrogram_metadata.eventList;
        numEvents = length(eventList);
        if isfield(region_power_spectrogram_metadata, 'Fs')
            Fs = region_power_spectrogram_metadata.Fs;
        else
            Fs = 496.031746;
        end
        if isfield(region_power_spectrogram_metadata, 'twin')
            twin = region_power_spectrogram_metadata.twin;
        else
            switch iTrialType
                case 1,
                    twin = [-1 2];
                otherwise,
                    twin = [-1 1];
            end
        end
        f = region_power_spectrogram_metadata.f;
        t = linspace(twin(1),twin(2),size(mean_spect,4));
        
        xticks = [twin(1),0,twin(2)];
        yticks = 0:20:100;
        figProps.n = numEvents;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;

%         numChannels_per_region   = zeros(numSessions, numRegions);
%         PROBABLY DON'T NEED ABOVE LINE FOR SUMMARY SHEET - ALREADY IN
%         METADATA?

        PDFname = fullfile(subject_powerSpectDir, [implantID '_' trialType '_powerSpect_summary.pdf']);

        numRegionPlots = 0;
        numPages = 0;
        for iRegion = 1 : length(ROI_list)
            
            % match the regions we want to plot with the appropriate
            % mean_spect array elements
            curRegion = ROI_list{iRegion};
            regionIdx = strcmpi(curRegion, region_power_spectrogram_metadata.regionList);
            
            rowNum = rem(iRegion, regions_per_page);
            if rowNum == 1
                h_fig = zeros(2,1);
                h_axes = zeros(2,figProps.m,figProps.n);
                for iPlotType = 1 : 2
                    [h_fig(iPlotType), b] = createFigPanels5(figProps);
                    if figProps.n == 1    % get a dim-2 matrix if number of columns = 1
                        h_axes(iPlotType,:) = b';
                    else
                        h_axes(iPlotType,:,:) = b;
                    end
                end
                page_regionList = curRegion;
                page_sessions_per_region = num2str(region_power_spectrogram_metadata.sessions_per_region(regionIdx));
                numPages = numPages + 1;
            else
                page_regionList = [page_regionList ', ' curRegion];
                page_sessions_per_region = [page_sessions_per_region ', ' ...
                                            num2str(region_power_spectrogram_metadata.sessions_per_region(regionIdx))];
            end
            if rowNum == 0; rowNum = regions_per_page; end
            
            if any(regionIdx)
                for iEventType = 1 : numEvents
                    
                    for iPlotType = 1 : 2
                        if figProps.n == 1    % get a dim-2 matrix if number of columns = 1
                            axes(h_axes(iPlotType, rowNum));
                        else
                            axes(h_axes(iPlotType, rowNum, iEventType));
                        end
                        switch iPlotType
                            case 1,
                                toPlot = squeeze(mean_spect(regionIdx,iEventType, :, :));
                                colorLim = power_clim;
                            case 2,
                                toPlot = log(squeeze(mean_spect(regionIdx,iEventType, :, :)));
                                colorLim = log_clim;
                        end
                        imagesc(t, f, toPlot);
                        set(gca,'ydir','normal');
                        set(gca,'clim',colorLim);
                        if rowNum == 1
                            title(eventList{iEventType});
                        end

                        if rowNum < figProps.m
                            set(gca,'xticklabel',{});
                        else
                            set(gca,'xtick',xticks);
                            xlabel('time (s)')
                        end
                        if iEventType > 1
                            set(gca,'yticklabel',{});
                        else
                            set(gca,'ytick',yticks);
                            ylabel('frequency (Hz)');
                        end
                        
                    end    % for iPlotType

                end    % for iEventType
            end
                
            if rem(iRegion, figProps.m) == 0 || iRegion == length(ROI_list)                % ADD IN TEXT HEADER HERE
                
                for iPlotType = 1 : 2
                    h_figAxes = createFigAxes(h_fig(iPlotType));
                    axes(h_figAxes);

                    switch iPlotType
                        case 1,
                            textStr{1} = ['power spectrograms averaged across sessions for ' implantID];
                            colorLim = power_clim;
                        case 2,
                            textStr{1} = ['log power spectrograms averaged across sessions for ' implantID];
                            colorLim = log_clim;
                    end
                    textStr{2} = ['Trial type: ' trialType];
                    textStr{3} = page_regionList;
                    textStr{4} = ['sessions for each region: ' page_sessions_per_region];
                    textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                    if numPages == 1 && iPlotType == 1
                        export_fig(PDFname, '-pdf', '-q101', '-painters','-nocrop');
                    else
                        export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
                    end
                    close(h_fig(iPlotType));
                end     % for iPlotType
                numPages = numPages + 1;
            end
                
        end    % for iRegion = 1 : length(ROI_list)
        
    end    % for iTrialType...
    
end    % for i_chDB...
            
            
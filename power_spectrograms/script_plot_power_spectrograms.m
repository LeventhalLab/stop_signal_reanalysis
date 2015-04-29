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
% ROI_list = {'EEG','Str','GP','STN','SNr'};
% ROI_idx = NaN(numSubs, length(ROI_list));
% numROI = zeros(1, length(ROI_list));
% numValidRegions = length(ROI_list);

numTrialTypes = length(trialTypeList);

% eventtWin(1,:) = [-1 2];   % for analysis of all trials
% % analysisWin(1) = 3;
% % stepSize(1)    = 3;
% for iTrialType = 2 : numTrialTypes
%     eventtWin(iTrialType,:) = [-1 1];
% %     analysisWin(iTrialType) = 0.1;
% %     stepSize(iTrialType)    = 0.05;
% end
channels_per_page = 5;
z_clim = [-3 3];
mrl_clim = [0 1e-2];
% colorLim = [-4 -1];

colorLim = [0, 1e-2];


figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;
figProps.m = channels_per_page;    % number of rows
sideMargins = 2; botMargin = 2.54;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
                                          
for i_chDB = 4 : 4%length(chDB_list)
    
%     if i_chDB < 7
%         hilbert_025Hz_directory = '/Volumes/RecordingsLeventhal2/stop-sig_reanalysis BU/Hilbert transformed LFP 025 Hz bins';
%     else
%         hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
%     end

    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
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
    channels = eval( chDB_info.name );
    
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
    channels = excludeChannels(cp, channels);
    
    sessionList = getSessionsfromChannelDB( channels );
    regionList  = getSubclassesfromChannelDB( channels );
    numSessions = length( sessionList );
    numRegions  = length( regionList );
    
    if i_chDB == 4
        startTrialType = 4;
    else
        startTrialType = 1;
    end
    for iTrialType = startTrialType : length(trialTypeList)

%         if i_chDB > 1 && (iTrialType == 2 || iTrialType == 3)
%             continue;
%         end
        trialType = trialTypeList{iTrialType}
        

%         powerSpect_metadata.eventList = eventList{iTrialType};
%         powerSpect_metadata.trialType = trialType;
%         powerSpect_metadata.eventtWin = twin;
        
        for i_testCh = 1 : length(channels)
            sample_powerSpect_name = ['powerSpect_' trialType '_' channels{i_testCh}.name '.mat'];
            test_session = channels{i_testCh}.session;
            powerSpect_sessionDir = fullfile(subject_powerSpectDir, test_session);
            sample_powerSpect_name = fullfile(powerSpect_sessionDir, sample_powerSpect_name);
            if exist(sample_powerSpect_name, 'file'); break; end
        end
        load(sample_powerSpect_name);    % sample file just to get the array sizes to initialize
        Fs = powerSpect_metadata.Fs;
        numEvents = size(powerSpect, 1);
        numFreqs  = size(powerSpect, 2);
        numSamps  = size(powerSpect, 3);
        
        twin = powerSpect_metadata.eventtWin;
        t = linspace(twin(1), twin(2), numSamps);
        xticks = [twin(1),0,twin(2)];
        f = powerSpect_metadata.freqs;
        yticks = 20:20:100;%min(f):20:max(f);
        
        eventList = powerSpect_metadata.eventList;

        figProps.n = numEvents;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;

        mean_session_spectrogram = zeros(numSessions, numRegions, numEvents, numFreqs, numSamps);
        numChannels_per_region   = zeros(numSessions, numRegions);

        startSession = 1;
        for iSession = startSession : numSessions
            
            powerSpect_sessionDir = fullfile(subject_powerSpectDir, sessionList{iSession});
            if ~exist(powerSpect_sessionDir, 'dir')
                mkdir(powerSpect_sessionDir);
            end

            cp = initChanParams();
            cp.session = sessionList{iSession};
            
            if ~isempty(strfind(trialType, 'nogo'))
                cp.task = 4;
            elseif ~isempty(strfind(trialType, 'stop'))
                cp.task = 3;
            else
                cp.task = -1;
            end
        
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);

            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
            sessionChannels = excludeChannels(cp, sessionChannels);

            cp = initChanParams();
            cp.tetrode = {'e2', 'e3', 'e03','e03'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels);continue;end
        
            numCh = length(sessionChannels);
                       
            numPages = 0;
            numChPlots = 0;
            PDFname = fullfile(powerSpect_sessionDir, [sessionList{iSession} '_' trialType '_power_spectrograms.pdf']);
            for iCh = 1 : numCh
%                 disp(sprintf('trialtype: %s, session %s, %d of %d; channel %d of %d', ...
%                     trialType, sessionList{iSession}, iSession, numSessions, iCh, numCh))

                ch = sessionChannels{iCh};
                if any(ch.wire.markedGood) == 0; continue; end

                powerSpect_name = ['powerSpect_' trialType '_' sessionChannels{iCh}.name '.mat'];
                powerSpect_name = fullfile(powerSpect_sessionDir, powerSpect_name);
                
                if ~exist(powerSpect_name,'file'); continue; end
                load(powerSpect_name);
                regionIdx = find(strcmpi(ch.location.subclass, regionList));
                if isempty(regionIdx); 
                    errorStr = sprintf('No region found for %s', ch.name);
                    error('plot_powerSpect:noregion',errorStr);
                end
                
                if size(mean_session_spectrogram, 3) == 1
                    mean_session_spectrogram(iSession, regionIdx, 1, :, :) = squeeze(mean_session_spectrogram(iSession, regionIdx, :, :, :)) + squeeze(powerSpect(1:numEvents, :, :));
                else
                    mean_session_spectrogram(iSession, regionIdx, :, :, :) = squeeze(mean_session_spectrogram(iSession, regionIdx, :, :, :)) + powerSpect(1:numEvents, :, :);
                end
                numChannels_per_region(iSession, regionIdx) = numChannels_per_region(iSession, regionIdx) + 1;
                numChPlots = numChPlots + 1;
                rowNum = rem(numChPlots, channels_per_page);
                if rowNum == 1
                    [h_fig, h_axes] = createFigPanels5(figProps);
                    page_chList = ch.name;
                    page_locList = ch.location.subclass;
                    numPages = numPages + 1;
                else
                    page_chList = [page_chList ', ' ch.name];
                    page_locList = [page_locList ', ' ch.location.subclass];
                end
                if rowNum == 0; rowNum = channels_per_page; end
                
                for iEventType = 1 : numEvents
                    axes(h_axes(rowNum, iEventType));

                    toPlot = squeeze(powerSpect(iEventType, :, :));
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

                end    % for iEventType
                
                if rem(numChPlots, figProps.m) == 0 || iCh == numCh                 % ADD IN TEXT HEADER HERE
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);

                    textStr{1} = 'power spectrograms';
                    textStr{2} = ['Trial type: ' trialType];
                    textStr{3} = page_chList;
                    textStr{4} = page_locList;
                    textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                    if numPages == 1
                        export_fig(PDFname, '-pdf', '-q101', '-painters','-nocrop');
                    else
                        export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
                    end
                    close(h_fig);
                end
                
            end    % for iCh = 1 : numCh
            
            % average across all tetrodes for each region for a single session
            numPagesForRegions = 0;
            for iRegion = 1 : numRegions
                if numChannels_per_region(iSession, iRegion) > 0
                    mean_session_spectrogram(iSession, iRegion, :, :, :) = mean_session_spectrogram(iSession, iRegion, :, :, :) / numChannels_per_region(iSession, iRegion);
                end
                if rem(iRegion, figProps.m) == 1
                    numPagesForRegions = numPagesForRegions + 1;
                    [h_fig, h_axes] = createFigPanels5(figProps);
                    page_regionList = regionList{iRegion};
                    page_chPerRegion = num2str(numChannels_per_region(iSession, iRegion));
                else
                    page_regionList = [page_regionList ', ' regionList{iRegion}];
                    page_chPerRegion = [page_chPerRegion ', ' num2str(numChannels_per_region(iSession, iRegion))];
                end
                
                rowNum = rem(iRegion, channels_per_page);
                if rowNum == 0; rowNum = channels_per_page; end

                for iEventType = 1 : numEvents
                    axes(h_axes(rowNum, iEventType));  

                    toPlot = squeeze(mean_session_spectrogram(iSession, iRegion, iEventType, :, :));
                    imagesc(t,f,toPlot);
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

                end
                
                if rem(iRegion, figProps.m) == 0 || iRegion == numRegions                 % ADD IN TEXT HEADER HERE
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);

                    textStr{1} = [sessionList{iSession} ', power spectrograms, single session averages within regions'];
                    textStr{2} = ['Trial type: ' trialType];
                    textStr{3} = page_regionList;
                    textStr{4} = ['number of channels per region: ' num2str(page_chPerRegion)];
                    textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                    export_fig(PDFname, '-pdf', '-q101', '-painters', '-append', '-nocrop');
                    close(h_fig);
                end

            end    % end for iRegion = 1 : numRegions

        end    % for iSession = 1 : numSessions

  %%  
        % average across all relevant regions for all sessions
        region_power_spectrogram_metadata.implantID = implantID;
        region_power_spectrogram_metadata.regionList = regionList;
        region_power_spectrogram_metadata.f = f;
        region_power_spectrogram_metadata.eventList = eventList;
        region_power_spectrogram_metadata.trialType = trialType;
        region_power_spectrogram_metadata.Fs = Fs;
        region_power_spectrogram_metadata.twin = twin;
              
        regionSummaryMatName = [implantID '_' trialType '_powerSpect_across_sessions.mat'];
        regionSummaryMatName = fullfile(subject_powerSpectDir, regionSummaryMatName);
    
        mean_spect = zeros(numRegions, numEvents, numFreqs, numSamps);
        numValidSessions = zeros(1, numRegions);
        numPagesForRegions = 0;
        PDFname = fullfile(subject_powerSpectDir, [implantID '_' trialType '_powerSpect_across_sessions.pdf']);
        for iRegion = 1 : numRegions

            for iSession = 1 : numSessions

                % mean_spect now contains the mean spectrogram matrix for
                % each region within each session
                if numChannels_per_region(iSession, iRegion) > 0
                    numValidSessions(iRegion) = numValidSessions(iRegion) + 1;
                    mean_spect(iRegion, :, :, :) = squeeze(mean_spect(iRegion, :, :, :)) + ...
                        squeeze(mean_session_spectrogram(iSession, iRegion, :, :, :));
                end
            end
            mean_spect(iRegion, :, :, :) = mean_spect(iRegion, :, :, :) / numValidSessions(iRegion);
    
            rowNum = rem(iRegion, channels_per_page);
            if rowNum == 0; rowNum = channels_per_page; end

            if rowNum == 1
                numPagesForRegions = numPagesForRegions + 1;
                [h_fig, h_axes] = createFigPanels5(figProps);
                page_regionList = regionList{iRegion};
                page_chPerRegion = num2str(numValidSessions(iRegion));
            else
                page_regionList = [page_regionList ', ' regionList{iRegion}];
                page_chPerRegion = [page_chPerRegion ', ' num2str(numValidSessions(iRegion))];
            end
        
            for iEventType = 1 : numEvents
                axes(h_axes(rowNum, iEventType));

                toPlot = squeeze(mean_spect(iRegion, iEventType, :, :));
                imagesc(t, f, toPlot);
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
                    set(gca,'ytick',yticks);
                end

            end    % for iEventType
        
            if rem(iRegion, figProps.m) == 0 || iRegion == numRegions
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);

                textStr{1} = [implantID ', power spectrograms, averages across sessions'];
                textStr{2} = ['Trial type: ' trialType];
                textStr{3} = page_regionList;
                textStr{4} = ['number of sessions per region: ' num2str(page_chPerRegion)];
                textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                if numPagesForRegions == 1
                    export_fig(PDFname, '-pdf', '-q101', '-painters','-nocrop');
                else
                    export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
                end
                close(h_fig);
            end

        end    % for iRegion = 1 : numRegions
        
        region_power_spectrogram_metadata.sessions_per_region = numValidSessions;
        save(regionSummaryMatName, 'mean_spect', 'region_power_spectrogram_metadata');
        
    end    % for iTrialType...
    
end    % for i_chDB...
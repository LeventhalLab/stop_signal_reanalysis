% script_plotRT_phaseCorrelations_gabors

% script to plot heat maps of the circular correlation coefficient between
% phase and RT

bitOrder = 'b';
colorLim = [0, 1];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
% phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_circstat';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabors';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabor_plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_phaseRT_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabors/IM166_phaseRTcorr_gabors/D2020091114/phase_RT_analysis_gabor_D2020091114T09.mat';

load(sample_phaseRT_file);
eventList = phaseRTcorr_metadata.eventList;
numEvents = length(eventList);
if numEvents == 6    % working around having included the same event twice in the phase-RT calculations
    if strcmpi(eventList{5},eventList{6})
        eventList = eventList(1:5);
        numEvents = 5;
    end
end
numSamps = size(circRTcorr, 3);
numFreqs = size(circRTcorr, 2);
twin = phaseRTcorr_metadata.twin;
t = linspace(twin(1), twin(2), numSamps);
f = phaseRTcorr_metadata.freqs;

xticks = [twin(1),0,twin(2)];
% yticks = 20:20:100;%min(f):20:max(f);
    
% try making a separate pdf for each session
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

desired_freq_ticks = [4,32,64,128,256,500];
                                          
for i_chDB = 1 : 1%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_phaseRTcorr_directory = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr_gabors']);
    if ~exist(subject_phaseRTcorr_directory, 'dir')
        disp([subject_phaseRTcorr_directory ' not found.']);
        continue;
    end

    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_gabor_RTcorr_plots']);
    if ~exist(subject_RTcorr_plots_directory, 'dir')
        mkdir(subject_RTcorr_plots_directory);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
    channels = excludeChannels(cp, channels);
    
    sessionList = getSessionsfromChannelDB( channels );
    regionList  = getSubclassesfromChannelDB( channels );
    
    numSessions = length( sessionList );
    numRegions  = length( regionList );
    
    mean_session_phaseRTcorr = zeros(numSessions, numRegions, numEvents, numFreqs, numSamps);
    numChannels_per_region = zeros(numSessions, numRegions);
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        numCh = length(sessionChannels);
        
        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorr_directory, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            disp([phaseRTcorr_sessionDir ' not found.']);
            continue;
        end

        numPages = 0;
        
        
        % find how many different regions there are for this set of
        % channels
        
        numChPlots = 0;
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            phaseRT_name = ['phase_RT_analysis_gabor_' ch.name '.mat'];
            phaseRT_name = fullfile(phaseRTcorr_sessionDir, phaseRT_name);
            
            if exist(phaseRT_name, 'file')
                load(phaseRT_name);
                % the variable RTphases is a cell array {iEvent, iFreq,
                % iQuant}, where iEvent is the event index, iFreq is the
                % frequency of interest index, and iQuant is the RT
                % quantile
                regionIdx = find(strcmpi(ch.location.subclass, regionList));
            else
                continue;
            end

            if isempty(regionIdx); 
                errorStr = sprintf('No region found for %s', ch.name);
                error('plotRTcorr:noregion',errorStr);
            end
            
            t_ticks = [phaseRTcorr_metadata.twin(1),0,phaseRTcorr_metadata.twin(2)];
            f = phaseRTcorr_metadata.freqs;
            freqTick_idx = zeros(1,length(desired_freq_ticks));
            for i_freqTick = 1 : length(desired_freq_ticks)
                freqTick_idx(i_freqTick) = find(abs(f - desired_freq_ticks(i_freqTick)) == ...  
                                                    min(abs(f - desired_freq_ticks(i_freqTick))));
            end
            
            
            mean_session_phaseRTcorr(iSession, regionIdx, :, :, :) = squeeze(mean_session_phaseRTcorr(iSession, regionIdx, :, :, :)) + circRTcorr(1:numEvents, :, :);
            numChannels_per_region(iSession, regionIdx) = numChannels_per_region(iSession, regionIdx) + 1;
            
            numChPlots = numChPlots + 1;
            rowNum = rem(numChPlots, channels_per_page);
            if rowNum == 1
                [h_fig, h_axes] = createFigPanels5(figProps);
                page_chList = ch.name;
                page_locList = ch.location.subclass;
                numPages = numPages + 1;
                pageName = sprintf('%s_phaseRTcorr_gabor_%02d.pdf',sessionList{iSession}, numPages);
%                 PDFname = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession} '_phase_RTcorr_gabor.pdf']);
                PDFname = fullfile(subject_RTcorr_plots_directory, pageName);
            else
                page_chList = [page_chList ', ' ch.name];
                page_locList = [page_locList ', ' ch.location.subclass];
            end
            if rowNum == 0; rowNum = channels_per_page; end
                
            for iEventType = 1 : numEvents
                axes(h_axes(rowNum, iEventType));
            
                x_ticks = t_ticks;
                y_ticks = freqTick_idx;
                toPlot = squeeze(circRTcorr(iEventType, :, :));
                imagesc(t, f, toPlot);
                set(gca,'ydir','normal',...
                        'clim',colorLim,...
                        'xtick',x_ticks,...
                        'ytick',round(f(y_ticks)));
%                 set(gca,'ydir','normal');
%                 set(gca,'clim',colorLim);
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
%                     set(gca,'ytick',yticks);
                    ylabel('frequency (Hz)');
                end
                
            end    % for iEventType
            
            if rem(numChPlots, figProps.m) == 0 || iCh == numCh                 % ADD IN TEXT HEADER HERE
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);
                
                textStr{1} = 'phase-RT correlations (circstat)';
                textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
                textStr{3} = page_chList;
                textStr{4} = page_locList;
                textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                                    
%                 if numPages == 1
%                     export_fig(PDFname, '-pdf', '-q101', '-painters');
%                 else
%                     export_fig(PDFname, '-pdf', '-q101', '-painters', '-append');
%                 end
                print(PDFname, '-dpdf');
                close(h_fig);
            end
            
        end    % for iCh = 1 : numCh

        % average across all tetrodes for each region for a single session
        numPagesForRegions = 0;
        for iRegion = 1 : numRegions
            if numChannels_per_region(iSession, iRegion) > 0
                mean_session_phaseRTcorr(iSession, iRegion, :, :, :) = mean_session_phaseRTcorr(iSession, iRegion, :, :, :) / numChannels_per_region(iSession, iRegion);
            end
            if rem(iRegion, figProps.m) == 1
                numPagesForRegions = numPagesForRegions + 1;
                [h_fig, h_axes] = createFigPanels5(figProps);
                page_regionList = regionList{iRegion};
                page_chPerRegion = num2str(numChannels_per_region(iSession, iRegion));
                pageName = sprintf('%s_phaseRTcorr_gabor_regionMean_%02d.pdf',sessionList{iSession},numPagesForRegions);
                PDFname = fullfile(subject_RTcorr_plots_directory, pageName);
            else
                page_regionList = [page_regionList ', ' regionList{iRegion}];
                page_chPerRegion = [page_chPerRegion ', ' num2str(numChannels_per_region(iSession, iRegion))];
            end
            
            rowNum = rem(iRegion, channels_per_page);
            if rowNum == 0; rowNum = channels_per_page; end
            
            for iEventType = 1 : numEvents
                axes(h_axes(rowNum, iEventType));  

                toPlot = squeeze(mean_session_phaseRTcorr(iSession, iRegion, iEventType, :, :));
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
                    set(gca,'ytick',round(f(y_ticks)));
                    ylabel('frequency (Hz)');
                end
                
            end
            
            if rem(iRegion, figProps.m) == 0 || iRegion == numRegions                 % ADD IN TEXT HEADER HERE
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);
                
                textStr{1} = [sessionList{iSession} ', phase-RT correlations (gabors circstat), single session averages within regions'];
                textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
                textStr{3} = page_regionList;
                textStr{4} = ['number of channels per region: ' num2str(page_chPerRegion)];
                textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                
                print(PDFname, '-dpdf');
                close(h_fig);
            end
            
        end    % end for iRegion = 1 : numRegions
        
    end    % for iSession = 1 : numSessions
  %%  
    % average across all regions for all sessions
    region_phase_RTcorr_metadata.implantID = implantID;
%     region_phase_RTcorr_metadata.regionList = regionList;
    region_phase_RTcorr_metadata.f = f;
    region_phase_RTcorr_metadata.t = t;
    region_phase_RTcorr_metadata.twin = twin;
    region_phase_RTcorr_metadata.eventList = eventList;
    
    regionSummaryMatName = [implantID '_phaseRTcorr_gabors_across_sessions.mat'];
    regionSummaryMatName = fullfile(subject_phaseRTcorr_directory, regionSummaryMatName);
    
    % don't bother writing zeros to disk for regions where there are no
    % sessions for this rat
    validRegions = false(1,numRegions);
    validRegionList = cell{1,1};
    numValidRegions = 0;
    for iRegion = 1 : numRegions
        if sum(squeeze(numChannels_per_region(:,iRegion))) > 0
            validRegions(iRegion) = true;
            numValidRegions = numValidRegions + 1;
            validRegionList{numValidRegions} = regionList{iRegion};
        end
    end
    region_phase_RTcorr_metadata.regionList = validRegionList;
    
    mean_phaseRTcorr = zeros(numValidRegions, numEvents, numFreqs, numSamps);
    numValidSessions = zeros(1, numValidRegions);
    numPagesForRegions = 0;
    PDFname = fullfile(subject_RTcorr_plots_directory, [implantID '_phaseRTcorr_gabors.pdf']);
    
    i_validRegion = 0;
    for iRegion = 1 : numRegions
        if rem(iRegion, figProps.m) == 1
            numPagesForRegions = numPagesForRegions + 1;
            [h_fig, h_axes] = createFigPanels5(figProps);
        end

        if ~validRegions(iRegion); continue; end
        i_validRegion = i_validRegion + 1;
        for iSession = 1 : numSessions
            
            % mean_session_phaseRTcorr now contains the mean phase-RT correlation matrix for
            % each region within each session
            if numChannels_per_region(iSession, iRegion) > 0
                numValidSessions(i_validRegion) = numValidSessions(i_validRegion) + 1;
                mean_phaseRTcorr(i_validRegion, :, :, :) = squeeze(mean_phaseRTcorr(i_validRegion, :, :, :)) + ...
                    squeeze(mean_session_phaseRTcorr(iSession, iRegion, :, :, :));
            end
        end
        mean_phaseRTcorr(i_validRegion, :, :, :) = mean_phaseRTcorr(i_validRegion, :, :, :) / numValidSessions(i_validRegion);
        
        rowNum = rem(i_validRegion, channels_per_page);
        if rowNum == 0; rowNum = channels_per_page; end

        if rowNum == 1
            numPagesForRegions = numPagesForRegions + 1;
            [h_fig, h_axes] = createFigPanels5(figProps);
            page_regionList = validRegionList{i_validRegion};
            page_chPerRegion = num2str(numValidSessions(i_validRegion));
            pageName = sprintf('%s_phaseRTcorr_gabor_regionSessionMean_%02d.pdf',sessionList{iSession},numPagesForRegions);
            PDFname = fullfile(subject_RTcorr_plots_directory, pageName);
        else
            page_regionList = [page_regionList ', ' validRegionList{i_validRegion}];
            page_chPerRegion = [page_chPerRegion ', ' num2str(numValidSessions(i_validRegion))];
        end
        
        for iEventType = 1 : numEvents
            axes(h_axes(rowNum, iEventType));
        
            toPlot = squeeze(mean_phaseRTcorr(i_validRegion, iEventType, :, :));
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
                set(gca,'ytick',round(f(y_ticks)));
            end
            
        end    % for iEventType

        if rem(i_validRegion, figProps.m) == 0 || iRegion == numRegions
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
        
    end    % for iRegion = 1 : numRegions
    region_phase_RTcorr_metadata.sessions_per_region = numValidSessions;
    save(regionSummaryMatName, 'mean_phaseRTcorr', 'region_phase_RTcorr_metadata');
    % 5/26/2016  NOTE, CHECK THAT mean_phaseRTcorr region indexing is
    % correct and that metadata are saved properly (sessions per region,
    % etc.)
end
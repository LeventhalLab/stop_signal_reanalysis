% script_plotRT_phaseCorrelations_20140506

% script to plot heat maps of the circular correlation coefficient between
% phase and RT

bitOrder = 'b';
yLim = [-1 1];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_circstat';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT correlation plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_phaseRT_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_circstat/IM166_phaseRTcorr_circstat/D2020091114/phase_RT_analysis_circstat_D2020091114T14.mat';

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
yticks = -1:1:1;%min(f):20:max(f);
    
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

targetFreq = 2.5;

for i_chDB = 1 : 4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_phaseRTcorr_directory = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr_circstat']);
    if ~exist(subject_phaseRTcorr_directory, 'dir')
        disp([subject_phaseRTcorr_directory ' not found.']);
        continue;
    end

    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_plots']);
    if ~exist(subject_RTcorr_plots_directory, 'dir')
        mkdir(subject_RTcorr_plots_directory);
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
    
    mean_session_phaseRTcorr = zeros(numSessions, numRegions, numEvents, numSamps);
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
        freqStr = [num2str(round(targetFreq)) 'Hz'];
        PDFname = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession} '_phase_RTcorr_' freqStr '.pdf']);
        
        % find how many different regions there are for this set of
        % channels
        
        numChPlots = 0;
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            phaseRT_name = ['phase_RT_analysis_circstat_' ch.name '.mat'];
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
            
            freqIdx = find(abs(phaseRTcorr_metadata.freqs - targetFreq) == min(abs(phaseRTcorr_metadata.freqs - targetFreq)));
            
            mean_session_phaseRTcorr(iSession, regionIdx, :, :) = squeeze(mean_session_phaseRTcorr(iSession, regionIdx, :, :)) + squeeze(circRTcorr(1:numEvents, freqIdx, :));
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
            
                toPlot = squeeze(circRTcorr(iEventType, freqIdx, :));
                plot(t, toPlot);
                hold on
                line([-1 1],[0 0],'color','k')
                line([0 0], [-1 1], 'color','k')
                set(gca,'ylim',[-1 1]);
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
                    ylabel('corr coef');
                end
                
            end    % for iEventType
            
            if rem(numChPlots, figProps.m) == 0 || iCh == numCh                 % ADD IN TEXT HEADER HERE
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);
                
                textStr{1} = 'phase-RT correlations (circstat)';
                textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
                textStr{3} = page_chList;
                textStr{4} = page_locList;
                textStr{5} = ['frequency: ' num2str(targetFreq) ' Hz'];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                                    
                if numPages == 1
                    export_fig(PDFname, '-pdf', '-q101', '-painters');
                else
                    export_fig(PDFname, '-pdf', '-q101', '-painters', '-append');
                end
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
            else
                page_regionList = [page_regionList ', ' regionList{iRegion}];
                page_chPerRegion = [page_chPerRegion ', ' num2str(numChannels_per_region(iSession, iRegion))];
            end
            
            rowNum = rem(iRegion, channels_per_page);
            if rowNum == 0; rowNum = channels_per_page; end
            
            for iEventType = 1 : numEvents
                axes(h_axes(rowNum, iEventType));  

                toPlot = squeeze(mean_session_phaseRTcorr(iSession, iRegion, iEventType, :));
                plot(t,toPlot);
                hold on
                line([-1 1],[0 0],'color','k')
                line([0 0], [-1 1], 'color','k')
                set(gca,'ylim',[-1 1]);
                
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
                    ylabel('corr coef');
                end
                
            end
            
            if rem(iRegion, figProps.m) == 0 || iRegion == numRegions                 % ADD IN TEXT HEADER HERE
                h_figAxes = createFigAxes(h_fig);
                axes(h_figAxes);
                
                textStr{1} = [sessionList{iSession} ', phase-RT correlations (circstat), single session averages within regions'];
                textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
                textStr{3} = page_regionList;
                textStr{4} = ['number of channels per region: ' num2str(page_chPerRegion)];
                textStr{5} = ['frequency: ' num2str(targetFreq) ' Hz'];
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                
                export_fig(PDFname, '-pdf', '-q101', '-painters', '-append');
                close(h_fig);
            end
            
        end    % end for iRegion = 1 : numRegions
        
    end    % for iSession = 1 : numSessions
  %%  
    % average across all regions for all sessions
    mean_phaseRTcorr = zeros(numRegions, numEvents, numSamps);
    numValidSessions = zeros(1, numRegions);
    numPagesForRegions = 0;
    PDFname = fullfile(subject_RTcorr_plots_directory, [implantID '_phaseRTcorr_' freqStr '.pdf']);
    for iRegion = 1 : numRegions
        if rem(iRegion, figProps.m) == 1
            numPagesForRegions = numPagesForRegions + 1;
            [h_fig, h_axes] = createFigPanels5(figProps);
        end

        for iSession = 1 : numSessions
            
            % mean_session_phaseRTcorr now contains the mean phase-RT correlation matrix for
            % each region within each session
            if numChannels_per_region(iSession, iRegion) > 0
                numValidSessions(iRegion) = numValidSessions(iRegion) + 1;
                mean_phaseRTcorr(iRegion, :, :) = squeeze(mean_phaseRTcorr(iRegion, :, :)) + ...
                    squeeze(mean_session_phaseRTcorr(iSession, iRegion, :, :));
            end
        end
        mean_phaseRTcorr(iRegion, :, :) = mean_phaseRTcorr(iRegion, :, :) / numValidSessions(iRegion);
        
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
        
            toPlot = squeeze(mean_phaseRTcorr(iRegion, iEventType, :));
            plot(t, toPlot);
            hold on
            line([-1 1],[0 0],'color','k')
            line([0 0], [-1 1], 'color','k')
            set(gca,'ylim',[-1 1]);
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
                ylabel('corr coef')
                set(gca,'ytick',yticks);
            end
            
        end    % for iEventType

        if rem(iRegion, figProps.m) == 0 || iRegion == numRegions
            h_figAxes = createFigAxes(h_fig);
            axes(h_figAxes);

            textStr{1} = [implantID ', phase-RT correlations (circstat), averages across sessions'];
            textStr{2} = ['Trial type: ' phaseRTcorr_metadata.trialType];
            textStr{3} = page_regionList;
            textStr{4} = ['number of sessions per region: ' num2str(page_chPerRegion)];
            textStr{5} = ['frequency: ' num2str(targetFreq) ' Hz'];
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

            if numPagesForRegions == 1
                export_fig(PDFname, '-pdf', '-q101', '-painters');
            else
                export_fig(PDFname, '-pdf', '-q101', '-painters', '-append');
            end
            close(h_fig);
        end
        
    end    % for iRegion = 1 : numRegions
    
end
% script_plotRT_powerCorrelations_gabor

% script to plot the correlation coefficients between power at different
% frequencies and RT

bitOrder = 'b';
colorLim = [0, 2e-5];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations_gabors';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT correlation gabor plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_powerRT_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations/IM164_powerRTcorr_gabor/D2220091006/power_RTcorr_D2220091006E01.mat';

load(sample_powerRT_file);
eventList = powerRTcorr_metadata.eventList;
numEvents = length(eventList);

% try making a separate pdf for each session
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page; figProps.n = numEvents;   % number of rows and columns, respectively

sideMargins = 1; topbotMargins = 2.54;

figProps.colSpacing = ones(1, figProps.n) * 0.5;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                             2 * sideMargins - ...
                                             sum(figProps.colSpacing)) / figProps.n;
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
    
    subject_powerRTcorr_directory = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr']);
    if ~exist(subject_powerRTcorr_directory, 'dir')
        disp([subject_powerRTcorr_directory ' not found.']);
        continue;
    end

    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_plots']);
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
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        numCh = length(sessionChannels);
        
        powerRTcorr_sessionDir = fullfile(subject_powerRTcorr_directory, sessionList{iSession});
        if ~exist(powerRTcorr_sessionDir, 'dir')
            disp([powerRTcorr_sessionDir ' not found.']);
            continue;
        end

        session_RTcorr_plots_directory = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession} '_RTcorr_gabor_plots']);
        if ~exist(session_RTcorr_plots_directory, 'dir')
            mkdir(session_RTcorr_plots_directory);
        end
    
        numPagesForSession = 0;
        PDFname = [sessionList{iSession} '_power_RTcorr_gabor'];
        
        % find how many different regions there are for this set of
        % channels
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            powerRT_name = ['power_RTcorr_' ch.name '.mat'];
            powerRT_name = fullfile(powerRTcorr_sessionDir, powerRT_name);
            
            if exist(powerRT_name, 'file')
                load(powerRT_name);
                % the variable RTphases is a cell array {iEvent, iFreq,
                % iQuant}, where iEvent is the event index, iFreq is the
                % frequency of interest index, and iQuant is the RT
                % quantile
                regionIdx = find(strcmpi(ch.location.subclass, regionList));
                
                rowIdx = (iCh - (numPagesForSession-1) * figProps.m);
                
            end
            
            mean_MIs(iSession, regionIdx, :, :) = squeeze(mean_MIs(iSession, regionIdx, :, :)) + MI;
            num_validMIs(iSession, regionIdx)   = num_validMIs(iSession, regionIdx) + 1;

            
            axes(h_axes(rowIdx, colIdx));
            
            imagesc(low_freqs, high_freqs, MI');
            set(gca,'ydir','normal');
            set(gca,'clim',colorLim);
            title([ch.name ', ' ch.location.name]);
            
            if rowIdx < figProps.m
                set(gca,'xticklabel',{});
            end
            if colIdx > 1
                set(gca,'yticklabel',{});
            end
            
            if rem(iCh, figProps.n * figProps.m) == 0 || iCh == numCh
                cur_PDFname = sprintf('%s_%02d.pdf',PDFname,numPagesForSession);
                cur_PDFname = fullfile(session_RTcorr_plots_directory, cur_PDFname);
                if numPagesForSession == 1
                    
                    % WORKING HERE - CHANGE TO PRINT
                    export_fig(pdfName, '-pdf', '-q101', '-painters');
                else
                    export_fig(pdfName, '-pdf', '-q101', '-painters', '-append');
                end
                close(h_fig);
            end
            
        end    % for iCh = 1 : numCh
        
        % average across all tetrodes for each region for a single session
        numPagesForRegions = 0;
        for iRegion = 1 : numRegions
            if num_validMIs(iSession, iRegion) > 0
                mean_MIs(iSession, iRegion, :, :) = mean_MIs(iSession, iRegion, :, :) / num_validMIs(iSession, iRegion);
            end
            if rem(iRegion, figProps.n * figProps.m) == 1
                numPagesForRegions = numPagesForRegions + 1;
                [h_fig, h_axes] = createFigPanels5(figProps);
            end
            
            rowIdx = ceil((iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) / figProps.n);
            colIdx = (iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) - ((rowIdx-1) * figProps.n);
            
            axes(h_axes(rowIdx, colIdx));
            
            imageMI = squeeze(mean_MIs(iSession, iRegion, :, :));
            imagesc(low_freqs, high_freqs, imageMI');
            set(gca,'ydir','normal');
            set(gca,'clim',colorLim);
            title([sessionList{iSession} ', ' regionList{iRegion} ', n = ' num2str(num_validMIs(iSession, iRegion))]);
            
            if rowIdx < figProps.m
                set(gca,'xticklabel',{});
            end
            if colIdx > 1
                set(gca,'yticklabel',{});
            end
            
            if rem(iRegion, figProps.n * figProps.m) == 0 || iRegion == numRegions
                export_fig(pdfName, '-pdf', '-q101', '-painters', '-append');
                close(h_fig);
            end
            
        end    % end for iRegion = 1 : numRegions
        
    end    % for iSession = 1 : numSessions
    
    % average across all regions for all sessions
    mean_regionMI = zeros(numRegions, size(MI,1), size(MI,2));
    numValidSessions = zeros(1, numRegions);
    numPagesForRegions = 0;
    pdfName = fullfile(subject_comodPlotsDir, [implantID '_phase_amp.pdf']);
    for iRegion = 1 : numRegions
        if rem(iRegion, figProps.n * figProps.m) == 1
            numPagesForRegions = numPagesForRegions + 1;
            [h_fig, h_axes] = createFigPanels5(figProps);
        end

        for iSession = 1 : numSessions
            
            % mean_MIs now contains the mean modulation index matrix for
            % each region within each session
            if num_validMIs(iSession, iRegion) > 0
                numValidSessions(iRegion) = numValidSessions(iRegion) + 1;
                mean_regionMI(iRegion, :, :) = squeeze(mean_regionMI(iRegion, :, :)) + ...
                    squeeze(mean_MIs(iSession, iRegion, :, :));
            end
        end
        mean_regionMI(iRegion, :, :) = mean_regionMI(iRegion, :, :) / numValidSessions(iRegion);
        
        rowIdx = ceil((iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) / figProps.n);
        colIdx = (iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) - ((rowIdx-1) * figProps.n);

        axes(h_axes(rowIdx, colIdx));
        
        imageMI = squeeze(mean_regionMI(iRegion, :, :));
        imagesc(low_freqs, high_freqs, imageMI');
        set(gca,'ydir','normal');
            set(gca,'clim',colorLim);
        title([implantID ', ' regionList{iRegion} ', n = ' num2str(numValidSessions(iRegion))]);
        
        if rowIdx < figProps.m
            set(gca,'xticklabel',{});
        end
        if colIdx > 1
            set(gca,'yticklabel',{});
        end

        if rem(iRegion, figProps.n * figProps.m) == 0 || iRegion == numRegions
            if numPagesForRegions == 1
                export_fig(pdfName, '-pdf', '-q101', '-painters');
            else
                export_fig(pdfName, '-pdf', '-q101', '-painters', '-append');
            end
            close(h_fig);
        end
        
    end    % for iRegion = 1 : numRegions
    
end
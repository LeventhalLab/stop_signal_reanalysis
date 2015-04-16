% script_plotRT_phaseCorrelations

% script to plot correlations between the phase of ongoing oscillations and
% RT.

% phase-RT correlations adapted from VanRullen et al, Ongoing
% EEG phase..., Frontiers in Psychology, 2011; details in Drewes and
% Vanrullen, "This is the rhtyhm of your eyes: the phase of ongoing
% electroencephalogram oscillations modulates saccadic reaction time.", J
% Neurosci, 2011
bitOrder = 'b';
colorLim = [0, 2e-5];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
phaseRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT correlation plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

% load a sample file
sample_phaseRT_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations/IM164_phaseRTcorr/D2220091005/phase_RT_analysis_D2220091005T02.mat';

load(sample_phaseRT_file);
eventList = phaseRTcorr_metadata.eventList;
numEvents = length(eventList);

% try making a separate pdf for each session
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = 2; figProps.n = numEvents;   % number of rows and columns, respectively

sideMargins = 1; topbotMargins = 2.54;

figProps.colSpacing = ones(1, figProps.n) * 0.5;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;

figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                             2 * sideMargins - ...
                                             sum(figProps.colSpacing)) / figProps.n;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              2 * topbotMargins - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_phaseRTcorr_directory = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr']);
    if ~exist(subject_phaseRTcorr_directory, 'dir')
        disp([subject_phaseRTcorr_directory ' not found.']);
        continue;
    end
    
    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_phaseRTcorr_plots']);
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

        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorr_directory, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            disp([phaseRTcorr_sessionDir ' not found.']);
            continue;
        end

        numPagesForSession = 0;
        PDFname = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession} '_phase_RTcorr.pdf']);
        
        % find how many different regions there are for this set of
        % channels
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            phaseRT_name = ['phase_RT_analysis_' ch.name '.mat'];
            phaseRT_name = fullfile(phaseRTcorr_sessionDir, phaseRT_name);
            
            if exist(phaseRT_name, 'file')
                if rem(iCh, figProps.m) == 1
                    numPagesForSession = numPagesForSession + 1;
                    [h_fig, h_axes] = createFigPanels2(figProps);
                end
                load(phaseRT_name);
                regionIdx = find(strcmpi(ch.location.subclass, regionList));
                
                rowIdx = (iCh - (numPagesForPhaseSession-1) * figProps.m);
                
                for iEvent = 1 : numEvents
                    colIdx = iEvent;
                    
                    
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
                if numPagesForSession == 1
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
                [h_fig, h_axes] = createFigPanels2(figProps);
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
            [h_fig, h_axes] = createFigPanels2(figProps);
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
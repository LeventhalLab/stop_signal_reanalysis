% script_plotRT_phaseCorrelations_gabor_20170105

% script to plot the correlation coefficients between phase at different
% frequencies and RT

colorLim = [0, 1];
desired_freq_ticks = [4,16,64,128,500];

chDB_directory         = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
phaseRTcorr_directory  = '/Volumes/Tbolt_02/stop-signal reanalysis/phase_RT_correlations_gabors';
RTcorr_plots_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/phase_RT correlation gabor plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_phaseRT_file = '/Volumes/Tbolt_02/stop-signal reanalysis/phase_RT_correlations_gabors/IM164_phaseRTcorr_gabors/D2220091006/phase_RT_analysis_gabor_D2220091006E01.mat';

load(sample_phaseRT_file);
eventList = phaseRTcorr_metadata.eventList;
numEvents = length(eventList);
f = phaseRTcorr_metadata.freqs;
f_to_plot = f;   % can change this later to be within a certain limit
numFreqs = length(f);
t = phaseRTcorr_metadata.t;
numSamps = length(t);
frange = [0,510];
desired_plot_f_ticks = desired_freq_ticks(desired_freq_ticks > frange(1) & ...
                                      desired_freq_ticks < frange(2));
freqTick_idx = zeros(1,length(desired_plot_f_ticks));
for i_freqTick = 1 : length(desired_plot_f_ticks)
    freqTick_idx(i_freqTick) = find(abs(f_to_plot - desired_plot_f_ticks(i_freqTick)) == ...  
                                        min(abs(f_to_plot - desired_plot_f_ticks(i_freqTick))));
end

% try making a separate pdf for each session
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page; figProps.n = numEvents;   % number of rows and columns, respectively

sideMargins = 1; botMargin = 2.54;

figProps.colSpacing = ones(1, figProps.n) * 0.5;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                             2 * sideMargins - ...
                                             sum(figProps.colSpacing)) / figProps.n;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

session_phaseRTcorr_metadata = phaseRTcorr_metadata;
for i_chDB = 1 : 4 %length(chDB_list)
    
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

    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_gabor_plots']);
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
        if isempty(sessionChannels);continue;end
        
        sessionRegionList  = getSubclassesfromChannelDB( sessionChannels );
        numSessionRegions = length(sessionRegionList);
        
        session_phaseRTcorr_metadata.regionList = sessionRegionList;
        
        numCh = length(sessionChannels);
        
        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorr_directory, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            disp([phaseRTcorr_sessionDir ' not found.']);
            continue;
        end

        session_RTcorr_plots_directory = fullfile(subject_RTcorr_plots_directory, [sessionList{iSession} '_RTcorr_gabor_plots']);
        if ~exist(session_RTcorr_plots_directory, 'dir')
            mkdir(session_RTcorr_plots_directory);
        end
    
        numPagesForSession = 0;
        PDFname = [sessionList{iSession} '_phase_RTcorr_gabor'];
        
        % find how many different regions there are for this set of
        % channels
        
        mean_phaseRTcorr = zeros(numSessionRegions, numEvents, numFreqs, numSamps);
        num_validCorrs = zeros(numSessionRegions, 1);
        
        numChPlots = 0;
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            phaseRT_name = ['phase_RT_analysis_gabor_' ch.name '.mat'];
            phaseRT_name = fullfile(phaseRTcorr_sessionDir, phaseRT_name);
            
%             phaseRTsurr_name = ['phase_RTcorr_' ch.name '_surr.mat'];
%             phaseRTsurr_name = fullfile(phaseRTcorr_sessionDir, phaseRTsurr_name);
            
            if exist(phaseRT_name, 'file')% && exist(phaseRTsurr_name, 'file')
                load(phaseRT_name);
%                 load(phaseRTsurr_name);
            else
                fprintf('%s does not exist, skipping...', phaseRT_name);
                continue;
            end
            
            regionIdx = find(strcmpi(ch.location.subclass, sessionRegionList));
            numChPlots = numChPlots + 1;
            rowNum = rem(iCh, channels_per_page);
            if rowNum == 0;rowNum = channels_per_page;end
            
            if rowNum == 1
                [h_phaseCorr_fig, h_phaseCorr_axes] = createFigPanels5(figProps);
%                 [h_phaseCorr_surr_fig, h_phaseCorr_surr_axes] = createFigPanels5(figProps);
                page_ChList = ch.name;
                page_regionList = ch.location.subclass;
%                 page_numSessionList = num2str(numRegionSessions(iRegion));
                numPagesForSession = numPagesForSession + 1;
            else
                page_ChList = [page_ChList ', ' ch.name];
                page_regionList = [page_regionList ', ' ch.location.subclass];
%                 page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
            end
        
            mean_phaseRTcorr(regionIdx, :, :, :) = squeeze(mean_phaseRTcorr(regionIdx, :, :, :)) + circRTcorr;
            num_validCorrs(regionIdx) = num_validCorrs(regionIdx) + 1;

            for iEvent = 1 : numEvents
                axes(h_phaseCorr_axes(rowNum, iEvent));

                toPlot = squeeze(circRTcorr(iEvent,:,:));
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
            
            if rem(iCh, figProps.m) == 0 || iCh == numCh
                h_figAxes = createFigAxes(h_phaseCorr_fig);

                textStr{1} = ['Circular correlation coefficient between phase and RT'];
                textStr{2} = page_ChList;
                textStr{3} = page_regionList;
%                 textStr{4} = ['number of channels per region: ' page_numChList];
                textStr{4} = sprintf('color limits: %d to %d',colorLim(1),colorLim(2));

                cur_PDFname = sprintf('%s_%02d.pdf',PDFname,numPagesForSession);
                cur_PDFname = fullfile(session_RTcorr_plots_directory, cur_PDFname);
                cur_figName = sprintf('%s_%02d.fig',PDFname,numPagesForSession);
                cur_figName = fullfile(session_RTcorr_plots_directory, cur_figName);

                axes(h_figAxes);
                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                print(cur_PDFname, '-dpdf');
                savefig(h_phaseCorr_fig,cur_figName,'compact');

                close(h_phaseCorr_fig);
            end
            
        end    % for iCh = 1 : numCh
        
        session_phaseRTcorr_metadata.numChPerRegion = num_validCorrs;
        
        for iRegion = 1 : numSessionRegions
            if num_validCorrs(iRegion) > 0
                mean_phaseRTcorr(iRegion, :, :, :) = mean_phaseRTcorr(iRegion, :, :, :) / num_validCorrs(iRegion);
            end
        end
            
        mean_phaseRTcorr_name = [PDFname '.mat'];
        mean_phaseRTcorr_name = fullfile(session_RTcorr_plots_directory, mean_phaseRTcorr_name);

        save(mean_phaseRTcorr_name, 'mean_phaseRTcorr', 'session_phaseRTcorr_metadata');
        
    end
    
end
            
%         
%         % average across all tetrodes for each region for a single session
%         numPagesForRegions = 0;
%         for iRegion = 1 : numSessionRegions
%             if num_validCorrs(iRegion) > 0
%                 mean_phaseRTcorr(iSession, iRegion, :, :) = mean_phaseRTcorr(iSession, iRegion, :, :) / num_validCorrs(iRegion);
%             end
%             
%             rowNum = rem(iRegion, figProps.m);
%             if rowNum == 0;rowNum = figProps.m;end
%             
%             if rowNum == 1
%                 numPagesForRegions = numPagesForRegions + 1;
%                 [h_phaseRTcorr_sessionRegion_fig, h_phaseRTcorr_sessionRegion_axes] = createFigPanels5(figProps);
%             end
%             
% %             rowIdx = ceil((iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) / figProps.n);
% %             colIdx = (iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) - ((rowIdx-1) * figProps.n);
%             
%             axes(h_phaseRTcorr_sessionRegion_axes(rowNum, iEvent));
%             
%             imageMI = squeeze(mean_MIs(iSession, iRegion, :, :));
%             imagesc(low_freqs, high_freqs, imageMI');
%             set(gca,'ydir','normal');
%             set(gca,'clim',colorLim);
%             title([sessionList{iSession} ', ' regionList{iRegion} ', n = ' num2str(num_validMIs(iSession, iRegion))]);
%             
%             if rowIdx < figProps.m
%                 set(gca,'xticklabel',{});
%             end
%             if colIdx > 1
%                 set(gca,'yticklabel',{});
%             end
%             
%             if rem(iRegion, figProps.n * figProps.m) == 0 || iRegion == numRegions
%                 export_fig(pdfName, '-pdf', '-q101', '-painters', '-append');
%                 close(h_fig);
%             end
%             
%         end    % end for iRegion = 1 : numRegions
%         
%     end    % for iSession = 1 : numSessions
%     
%     % average across all regions for all sessions
%     mean_regionMI = zeros(numRegions, size(MI,1), size(MI,2));
%     numValidSessions = zeros(1, numRegions);
%     numPagesForRegions = 0;
%     pdfName = fullfile(subject_comodPlotsDir, [implantID '_phase_amp.pdf']);
%     for iRegion = 1 : numRegions
%         if rem(iRegion, figProps.n * figProps.m) == 1
%             numPagesForRegions = numPagesForRegions + 1;
%             [h_fig, h_axes] = createFigPanels5(figProps);
%         end
% 
%         for iSession = 1 : numSessions
%             
%             % mean_MIs now contains the mean modulation index matrix for
%             % each region within each session
%             if num_validMIs(iSession, iRegion) > 0
%                 numValidSessions(iRegion) = numValidSessions(iRegion) + 1;
%                 mean_regionMI(iRegion, :, :) = squeeze(mean_regionMI(iRegion, :, :)) + ...
%                     squeeze(mean_MIs(iSession, iRegion, :, :));
%             end
%         end
%         mean_regionMI(iRegion, :, :) = mean_regionMI(iRegion, :, :) / numValidSessions(iRegion);
%         
%         rowIdx = ceil((iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) / figProps.n);
%         colIdx = (iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) - ((rowIdx-1) * figProps.n);
% 
%         axes(h_axes(rowIdx, colIdx));
%         
%         imageMI = squeeze(mean_regionMI(iRegion, :, :));
%         imagesc(low_freqs, high_freqs, imageMI');
%         set(gca,'ydir','normal');
%             set(gca,'clim',colorLim);
%         title([implantID ', ' regionList{iRegion} ', n = ' num2str(numValidSessions(iRegion))]);
%         
%         if rowIdx < figProps.m
%             set(gca,'xticklabel',{});
%         end
%         if colIdx > 1
%             set(gca,'yticklabel',{});
%         end
% 
%         if rem(iRegion, figProps.n * figProps.m) == 0 || iRegion == numRegions
%             if numPagesForRegions == 1
%                 export_fig(pdfName, '-pdf', '-q101', '-painters');
%             else
%                 export_fig(pdfName, '-pdf', '-q101', '-painters', '-append');
%             end
%             close(h_fig);
%         end
%         
%     end    % for iRegion = 1 : numRegions
%     
% end
% script_plotRT_powerCorrelations_gabor_20170105

% script to plot the correlation coefficients between power at different
% frequencies and RT

colorLim = [-1, 1];
desired_freq_ticks = [4,16,64,128,500];

chDB_directory         = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
powerRTcorr_directory  = '/Volumes/Tbolt_02/stop-signal reanalysis/power_RT_correlations_gabors';
RTcorr_plots_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/power_RT correlation gabor plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_powerRT_file = '/Volumes/Tbolt_02/stop-signal reanalysis/power_RT_correlations_gabors/IM164_powerRTcorr_gabor/D2220091006/power_RTcorr_D2220091006E01.mat';

load(sample_powerRT_file);
eventList = powerRTcorr_metadata.eventList;
numEvents = length(eventList);
f = powerRTcorr_metadata.freqs;
f_to_plot = f;   % can change this later to be within a certain limit
numFreqs = length(f);
t = powerRTcorr_metadata.t;
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

session_powerRTcorr_metadata = powerRTcorr_metadata;
for i_chDB = 1 : 4 %length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_powerRTcorr_directory = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr_gabor']);
    if ~exist(subject_powerRTcorr_directory, 'dir')
        disp([subject_powerRTcorr_directory ' not found.']);
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
        
        session_powerRTcorr_metadata.regionList = sessionRegionList;
        
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
        
        mean_powerRTcorr = zeros(numSessionRegions, numEvents, numFreqs, numSamps);
        num_validCorrs = zeros(numSessionRegions, 1);
        
        numChPlots = 0;
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            powerRT_name = ['power_RTcorr_' ch.name '.mat'];
            powerRT_name = fullfile(powerRTcorr_sessionDir, powerRT_name);
            
%             powerRTsurr_name = ['power_RTcorr_' ch.name '_surr.mat'];
%             powerRTsurr_name = fullfile(powerRTcorr_sessionDir, powerRTsurr_name);
            
            if exist(powerRT_name, 'file')% && exist(powerRTsurr_name, 'file')
                load(powerRT_name);
%                 load(powerRTsurr_name);
            else
                fprintf('%s does not exist, skipping...\n', powerRT_name);
                if iCh < numCh || ~isvalid(h_powerCorr_fig)    % if the figure was already deleted, don't save and delete it again
                    continue;
                else    % need to close out the previous figure and save it
                    h_figAxes = createFigAxes(h_powerCorr_fig);

                    textStr{1} = ['Spearman correlation coefficient between power and RT'];
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
                    savefig(h_powerCorr_fig,cur_figName,'compact');

                    close(h_powerCorr_fig);
                    continue;
                end
            end
            
            regionIdx = find(strcmpi(ch.location.subclass, sessionRegionList));
            numChPlots = numChPlots + 1;
            rowNum = rem(numChPlots, channels_per_page);
            if rowNum == 0;rowNum = channels_per_page;end
            
            if rowNum == 1
                [h_powerCorr_fig, h_powerCorr_axes] = createFigPanels5(figProps);
%                 [h_powerCorr_surr_fig, h_powerCorr_surr_axes] = createFigPanels5(figProps);
                page_ChList = ch.name;
                page_regionList = ch.location.subclass;
%                 page_numSessionList = num2str(numRegionSessions(iRegion));
                numPagesForSession = numPagesForSession + 1;
            else
                page_ChList = [page_ChList ', ' ch.name];
                page_regionList = [page_regionList ', ' ch.location.subclass];
%                 page_numSessionList = [page_numSessionList ', ' num2str(numRegionSessions(iRegion))];
            end
        
            mean_powerRTcorr(regionIdx, :, :, :) = squeeze(mean_powerRTcorr(regionIdx, :, :, :)) + powerRTcorr;
            num_validCorrs(regionIdx) = num_validCorrs(regionIdx) + 1;

            for iEvent = 1 : numEvents
                axes(h_powerCorr_axes(rowNum, iEvent));

                toPlot = squeeze(powerRTcorr(iEvent,:,:));
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
            
            if rem(numChPlots, figProps.m) == 0 || iCh == numCh
                h_figAxes = createFigAxes(h_powerCorr_fig);

                textStr{1} = ['Spearman correlation coefficient between power and RT'];
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
                savefig(h_powerCorr_fig,cur_figName,'compact');

                close(h_powerCorr_fig);
            end
            
        end    % for iCh = 1 : numCh
        
        session_powerRTcorr_metadata.numChPerRegion = num_validCorrs;
        
        for iRegion = 1 : numSessionRegions
            if num_validCorrs(iRegion) > 0
                mean_powerRTcorr(iRegion, :, :, :) = mean_powerRTcorr(iRegion, :, :, :) / num_validCorrs(iRegion);
            end
        end
            
        mean_powerRTcorr_name = [PDFname '.mat'];
        mean_powerRTcorr_name = fullfile(session_RTcorr_plots_directory, mean_powerRTcorr_name);

        save(mean_powerRTcorr_name, 'mean_powerRTcorr', 'session_powerRTcorr_metadata');
        
    end
    
end
            
%         
%         % average across all tetrodes for each region for a single session
%         numPagesForRegions = 0;
%         for iRegion = 1 : numSessionRegions
%             if num_validCorrs(iRegion) > 0
%                 mean_powerRTcorr(iSession, iRegion, :, :) = mean_powerRTcorr(iSession, iRegion, :, :) / num_validCorrs(iRegion);
%             end
%             
%             rowNum = rem(iRegion, figProps.m);
%             if rowNum == 0;rowNum = figProps.m;end
%             
%             if rowNum == 1
%                 numPagesForRegions = numPagesForRegions + 1;
%                 [h_powerRTcorr_sessionRegion_fig, h_powerRTcorr_sessionRegion_axes] = createFigPanels5(figProps);
%             end
%             
% %             rowIdx = ceil((iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) / figProps.n);
% %             colIdx = (iRegion - (numPagesForRegions-1) * figProps.m * figProps.n) - ((rowIdx-1) * figProps.n);
%             
%             axes(h_powerRTcorr_sessionRegion_axes(rowNum, iEvent));
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
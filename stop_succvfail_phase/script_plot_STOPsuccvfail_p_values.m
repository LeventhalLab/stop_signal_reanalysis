% script_plot_STOPsuccvfail_p_values

% script to load in comparisons between STOP-success and fail trials
% (phase) and make plots of the p-values

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

channels_per_page = 5;
colorLim = [0 1];
vecDiff_colorLim = [0 2];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
                                          
makePlots = true;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        disp([subject_stopPhaseDir ' not found. Skipping ' implantID '...'])
        continue
    end

    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.task = 3;    % only stop-signal sessions
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        regionList = getRegionsfromChannelDB(sessionChannels);
        numRegions = length(regionList);
        
        totalSessionChannels = length(sessionChannels);
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
        if ~exist(session_stopPhaseDir, 'dir')
            disp([session_stopPhaseDir ' not found. Skipping ' sessionList{iSession} '...'])
            continue
        end
        
        boot_region_mean_saveName = ['bootStrapped_p_stopSuccvFail_' sessionList{iSession} '.mat'];
        boot_region_mean_saveName = fullfile(session_stopPhaseDir, boot_region_mean_saveName);
        
        ww_region_mean_saveName = ['ww_p_stopSuccvFail_' sessionList{iSession} '.mat'];
        ww_region_mean_saveName = fullfile(session_stopPhaseDir, ww_region_mean_saveName);
       
        boot_fig_saveName = ['bootStrapped_p_stopSuccvFail_' sessionList{iSession} '.pdf'];
        boot_fig_saveName = fullfile(session_stopPhaseDir, boot_fig_saveName);
        
        ww_fig_saveName = ['ww_p_stopSuccvFail_' sessionList{iSession} '.pdf'];
        ww_fig_saveName = fullfile(session_stopPhaseDir, ww_fig_saveName);
        
        vecDiff_fig_saveName = ['vecDiff_stopSuccvFail_' sessionList{iSession} '.pdf'];
        vecDiff_fig_saveName = fullfile(session_stopPhaseDir, vecDiff_fig_saveName);
        
        for iCh = 1 : totalSessionChannels
            wwtest_fileName = ['STOPsucc_fail_wwtest_' sessionChannels{iCh}.name '.mat'];
            wwtest_fileName = fullfile(session_stopPhaseDir, wwtest_fileName);
            if ~exist(wwtest_fileName, 'file'); continue; end
            
            load(wwtest_fileName);
            break;
        end
        numEventTypes = length(STOPanal_metadata.eventList);
        numFreqs  = size(STOPanal_metadata.freqBands, 1);
        numSamps  = size(ww_p, 3);
        twin = STOPanal_metadata.twin;
        f = mean(STOPanal_metadata.freqBands, 2);
        
        t = linspace(twin(1), twin(2), numSamps);

        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
        numPages = 0;
        wwp_by_region = cell(1, numRegions);    %zeros(numRegions, numEventTypes, numFreqs, numSamps);
        bootp_by_region = cell(1, numRegions);    %zeros(numRegions, numEventTypes, numFreqs, numSamps);
        vecDiff_by_region = cell(1, numRegions);

        num_ch_per_region = zeros(1, numRegions);
        numChPlots = 0;
        for iRegion = 1 : numRegions
                
            cp = initChanParams();
            cp.locationName = regionList{iRegion};            
            chList = extractChannels(cp, sessionChannels);
            if isempty(chList); continue; end
            regionChannels = sessionChannels(chList);
            numCh = length(regionChannels);
            
            wwp_by_region{iRegion} = zeros(numCh, numEventTypes, numFreqs, numSamps);
            bootp_by_region{iRegion} = zeros(numCh, numEventTypes, numFreqs, numSamps);
            vecDiff_by_region{iRegion} = zeros(numCh, numEventTypes, numFreqs, numSamps);
            mean_wwp = zeros(numRegions, numEventTypes, numFreqs, numSamps);
            mean_bootp = zeros(numRegions, numEventTypes, numFreqs, numSamps);
            mean_vecDiff = zeros(numRegions, numEventTypes, numFreqs, numSamps);
            
            num_ch_per_region(iRegion) = numCh;
            
            for iCh = 1 : numCh
                iCh
            
                ch = regionChannels{iCh};
                wwtest_fileName = ['STOPsucc_fail_wwtest_' ch.name '.mat'];
                wwtest_fileName = fullfile(session_stopPhaseDir, wwtest_fileName);
                bootstrap_fileName = ['STOPsucc_fail_boot_' ch.name '.mat'];
                bootstrap_fileName = fullfile(session_stopPhaseDir, bootstrap_fileName);
                if ~exist(wwtest_fileName, 'file'); continue; end

                load(wwtest_fileName);
                load(bootstrap_fileName);
                
                numChPlots = numChPlots + 1;
                rowNum = rem(numChPlots, channels_per_page);
                if rowNum == 1
                    if makePlots
                        [h_ww_fig, h_ww_axes] = createFigPanels5(figProps);
                        [h_boot_fig, h_boot_axes] = createFigPanels5(figProps);
                        [h_vecDiff_fig, h_vecDiff_axes] = createFigPanels5(figProps);
                    end
                    page_chList = ch.name;
                    page_locList = ch.location.subclass;
                    numPages = numPages + 1;
                else
                    page_chList = [page_chList ', ' ch.name];
                    page_locList = [page_locList ', ' ch.location.subclass];
                end
                if rowNum == 0; rowNum = channels_per_page; end
                
                wwp_by_region{iRegion}(iCh, :, :, :) = ww_p;
                bootp_by_region{iRegion}(iCh, :, :, :) = boot_p;
                
                % calculate absolute distance between mean resultant
                % vectors foir stop-success and fail trials. r(1, :, :, ;)
                % is for stop-success
                r_succ = squeeze(r(1, :, :, :)); r_fail = squeeze(r(2, :, :, :));
                theta_succ = squeeze(theta(1, :, :, :)); theta_fail = squeeze(theta(2, :, :, :));
                vecDiff =  abs(r_succ .* exp(1i*theta_succ) - r_fail .* exp(1i*theta_fail));
                vecDiff_by_region{iRegion}(iCh, :, :, :) = vecDiff;
            
                if makePlots
                    for iEventType = 1 : numEventTypes
                    
                        ww_toPlot = squeeze(ww_p(iEventType, :, :))';
                        boot_toPlot = squeeze(boot_p(iEventType, :, :))';
                        vecDiff_toPlot = squeeze(vecDiff(iEventType, :, :))';

                        axes(h_ww_axes(rowNum, iEventType));
                        imagesc(t, f, ww_toPlot);
                        set(gca,'ydir','normal','clim',colorLim);
                        if rowNum == 1
                            title(STOPanal_metadata.eventList{iEventType});
                        end
                        if rowNum < channels_per_page
                            set(gca,'xticklabel',[]);
                        end
                        if iEventType > 1
                            set(gca,'yticklabel',[]);
                        end
                        
                        axes(h_boot_axes(rowNum, iEventType));
                        imagesc(t, f, boot_toPlot);
                        set(gca,'ydir','normal','clim',colorLim);
                        if rowNum == 1
                            title(STOPanal_metadata.eventList{iEventType});
                        end
                        if rowNum < channels_per_page
                            set(gca,'xticklabel',[]);
                        end
                        if iEventType > 1
                            set(gca,'yticklabel',[]);
                        end
                        
                        axes(h_vecDiff_axes(rowNum, iEventType));
                        imagesc(t, f, vecDiff_toPlot);
                        set(gca,'ydir','normal','clim',vecDiff_colorLim);
                        if rowNum == 1
                            title(STOPanal_metadata.eventList{iEventType});
                        end
                        if rowNum < channels_per_page
                            set(gca,'xticklabel',[]);
                        end
                        if iEventType > 1
                            set(gca,'yticklabel',[]);
                        end
                        
                    end    % for iEventType = 1 : numEventTypes
                        
                    if (rowNum == channels_per_page || numChPlots == totalSessionChannels)
                        h_figAxes = createFigAxes(h_ww_fig);
                        axes(h_figAxes);

                        textStr{1} = ['p-values for stop-success vs failure comparisons, Watson-Williams test'];
                        textStr{2} = page_chList;
                        textStr{3} = page_locList;
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                        if numPages == 1
                            export_fig(ww_fig_saveName,'-pdf','-q101','-painters');
                        else
                            export_fig(ww_fig_saveName,'-pdf','-q101','-painters','-append');
                        end
                        close(h_ww_fig);


                        h_figAxes = createFigAxes(h_boot_fig);
                        axes(h_figAxes);

                        textStr{1} = ['p-values for stop-success vs failure comparisons, bootstrapping test'];
                        textStr{2} = page_chList;
                        textStr{3} = page_locList;
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                        if numPages == 1
                            export_fig(boot_fig_saveName,'-pdf','-q101','-painters');
                        else
                            export_fig(boot_fig_saveName,'-pdf','-q101','-painters','-append');
                        end
                        close(h_boot_fig);
                        
                        h_figAxes = createFigAxes(h_vecDiff_fig);
                        axes(h_figAxes);

                        textStr{1} = ['absolute difference in mean resultant vectors'];
                        textStr{2} = page_chList;
                        textStr{3} = page_locList;
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                        if numPages == 1
                            export_fig(vecDiff_fig_saveName,'-pdf','-q101','-painters');
                        else
                            export_fig(vecDiff_fig_saveName,'-pdf','-q101','-painters','-append');
                        end
                        close(h_vecDiff_fig);

                    end
                    
                end    % if makePlots
                
            end    % for iCh...
            
            mean_wwp(iRegion, :, : ,:) = squeeze(mean(wwp_by_region{iRegion}, 1));
            mean_bootp(iRegion, :, : ,:) = squeeze(mean(bootp_by_region{iRegion}, 1));
            
        end    % for iRegion...
        STOPregion_metadata = STOPanal_metadata;
        STOPregion_metadata.regionList = regionList;
        STOPregion_metadata.num_ch_per_region = num_ch_per_region;
        STOPregion_metadata.f = f;
        save(ww_region_mean_saveName, 'mean_wwp', 'STOPregion_metadata');
        save(boot_region_mean_saveName, 'mean_bootp', 'STOPregion_metadata');
        
        % create a new page of figures to show the region means for
        % each session
        if makePlots
            for iRegion = 1 : numRegions
                rowNum = rem(iRegion, channels_per_page);
                if rowNum == 1
                    [h_fig, h_axes] = createFigPanels5(figProps);
                    page_regionList = regionList{iRegion};
                    page_numChList  = ['Number of Channels per Region: ' num2str(num_ch_per_region(iRegion))];
                else
                    page_regionList = [page_regionList ', ' regionList{iRegion}];
                    page_numChList  = [page_numChList ', ' num2str(num_ch_per_region(iRegion))];
                end
        
                    % WORKING HERE...

                for iFreq = STOPanal_metadata.freqComplete + 1 : numFreqs
                    iFreq
                    for i_t = 1 : numSamps
                        
                        phases_correct = squeeze(correctSTOP_phase(iEvent, iFreq, :, i_t));
                        phases_failed  = squeeze(failedSTOP_phase(iEvent, iFreq, :, i_t));
                        [ww_p(iEvent, iFreq, i_t), ~] = circ_wwtest(phases_correct, phases_failed);
                        temp = boot_angle_test(phases_correct, phases_failed, 'nboot', nBoot);
                        boot_p(iEvent, iFreq, i_t) = boot_angle_test(phases_correct, phases_failed, 'nboot', nBoot);
                        for i_outcome = 1 : 2
                            r(i_outcome, iEvent, iFreq, i_t) = abs(mean(exp(1i * phases_correct)));
                            theta(i_outcome, iEvent, iFreq, i_t) = angle(mean(exp(1i * phases_correct)));
                        end
                        
                    end
                    STOPanal_metadata.freqComplete = STOPanal_metadata.freqComplete + 1;
                    save(wwtest_saveName, 'ww_p', 'r', 'theta', 'STOPanal_metadata');
                    save(bootstrap_saveName, 'boot_p', 'r', 'theta', 'STOPanal_metadata');
                end
                STOPanal_metadata_loaded = false;
                toc
            end    % for iEvent...
            
        end
        
    end
    
end
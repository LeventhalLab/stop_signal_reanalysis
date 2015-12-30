% script_review_phaseAmp_Gabors

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed_Gabors';

makePlots = true;
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};
numTrialTypes = length(trialTypeList);

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
regions_per_page = 5;
z_clim = [-4 4];
mrl_clim = [0 1e-2];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = regions_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;

figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

var_to_plot = {'mrl','zscore'};
plotTypes = {'const_phase_f', 'const_amp_f', 'averaged_t'};
phase_freq = [2.0, 4.0, 8.0];       % Hz, for the time-frequency plots, can only look at one "phase" frequency at a time
amp_freq   = [19.5, 49.5];          % Hz, for the time-frequency plots, can only look at one "amplitude" frequency at a time
plotFreqs{1} = phase_freq;
plotFreqs{2} = amp_freq;
plotFreqs{3} = 1;    % dummy variable to make the indexing work later
totalPlotTypes = length(var_to_plot) * (1 + length(phase_freq) + length(amp_freq));
% numPlotTypes = length(plotTypes);

desired_amp_freq_ticks = [10,20,50,80];    % check that these frequencies are the ones actually being marked
desired_phase_freq_ticks = [2,4,8,16,20];

for i_chDB = 1 : 1%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        continue;
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);  
    
    for iTrialType = 1 : 1%length(trialTypeList)
        trialType = trialTypeList{iTrialType};

        fig_saveName = cell(length(var_to_plot), length(plotTypes));
        mat_saveName = [implantID '_phaseAmp_' trialTypeList{iTrialType} '_Gabor_summary.mat'];
        subject_trialType_dir = fullfile(subject_phaseAmpdir, [implantID '_phaseAmp_' trialType '_Gabors']);
        if ~exist(subject_trialType_dir,'dir')
            mkdir(subject_trialType_dir);
        end
        mat_saveName = fullfile(subject_trialType_dir, mat_saveName);
        if exist(mat_saveName,'file');continue;end
        
        mat_z_saveName = [implantID '_phaseAmp_' trialTypeList{iTrialType} '_z_Gabor_summary.mat'];
        mat_z_saveName = fullfile(subject_trialType_dir, mat_z_saveName);
        for iVar = 1 : length(var_to_plot)
            for iPlotType = 1 : length(plotTypes)
                fig_saveName{iVar, iPlotType} = cell(1,length(plotFreqs{iPlotType}));

                for iFreq = 1 : length(plotFreqs{iPlotType})
                    freqStr = sprintf('%04.0f',plotFreqs{iPlotType}(iFreq));
                    fig_saveName{iVar, iPlotType}{iFreq} = ...
                        ['phaseAmp_' trialTypeList{iTrialType} '_' var_to_plot{iVar} '_' plotTypes{iPlotType} '_' freqStr '_Gabor.pdf'];
                    fig_saveName{iVar, iPlotType}{iFreq} = fullfile(subject_trialType_dir, fig_saveName{iVar, iPlotType}{iFreq});
                end

            end    % for iPlotType...
        end
        
        % load the first session mean mrl file to figure out how to
        % allocate memory
        for iSession = 1 : numSessions
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            region_meanName = ['gaborPhaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_meanName = fullfile(phaseAmp_sessionDir, region_meanName);
            if ~exist(region_meanName, 'file'); continue; end
            load(region_meanName);
            break;
        end
        
        phaseAmpSummary_metadata.phase_f = region_phaseAmp_metadata.phase_f;
        phaseAmpSummary_metadata.amp_f = region_phaseAmp_metadata.amp_f;
        phaseAmpSummary_metadata.eventList = region_phaseAmp_metadata.eventList;
        phaseAmpSummary_metadata.eventtWin = region_phaseAmp_metadata.eventtWin;
        phaseAmpSummary_metadata.Fs = region_phaseAmp_metadata.Fs;
        phaseAmpSummary_metadata.regionList = region_phaseAmp_metadata.regionList;
        
        eventtWin = region_phaseAmp_metadata.eventtWin;
        numSamps = size(mean_mrl, 5);
        t = linspace(eventtWin(1), eventtWin(2), numSamps);
        
        phase_f = region_phaseAmp_metadata.phase_f;
        amp_f = region_phaseAmp_metadata.amp_f;
        
        numEventTypes = length(region_phaseAmp_metadata.eventList);
        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
                                                 
        phase_freq_idx = zeros(size(phase_freq));
        amp_freq_idx = zeros(size(amp_freq));
        for iFreq = 1 : length(phase_freq)
            phase_freq_idx(iFreq) = find(abs(phase_f - phase_freq(iFreq)) == ...
                                         min(abs(phase_f - phase_freq(iFreq))));
        end
        for iFreq = 1 : length(amp_freq)
            amp_freq_idx(iFreq) = find(abs(amp_f - amp_freq(iFreq)) == ...
                                       min(abs(amp_f - amp_freq(iFreq))));
        end
        
        amp_freqTick_idx = zeros(1, length(desired_amp_freq_ticks));
        phase_freqTick_idx = zeros(1, length(desired_phase_freq_ticks));
        amp_freqTick_label = zeros(1, length(desired_amp_freq_ticks));
        phase_freqTick_label = zeros(1, length(desired_phase_freq_ticks));
        for i_freqTick = 1 : length(desired_amp_freq_ticks)
            amp_freqTick_idx(i_freqTick) = find(abs(amp_f - desired_amp_freq_ticks(i_freqTick)) == ...
                                                min(abs(amp_f - desired_amp_freq_ticks(i_freqTick))));
            amp_freqTick_label(i_freqTick) = round(amp_f(amp_freqTick_idx(i_freqTick)));
        end
        for i_freqTick = 1 : length(desired_phase_freq_ticks)
            phase_freqTick_idx(i_freqTick) = find(abs(phase_f - desired_phase_freq_ticks(i_freqTick)) == ...
                                                min(abs(phase_f - desired_phase_freq_ticks(i_freqTick))));
            phase_freqTick_label(i_freqTick) = round(phase_f(phase_freqTick_idx(i_freqTick)));
        end
        t_ticks = [region_phaseAmp_metadata.eventtWin(1),0,region_phaseAmp_metadata.eventtWin(2)];
        
        % NEED TO FIGURE OUT HOW TO DO THESE CALCULATIONS PIECE-WISE TO
        % PRESERVE MEMORY; WILL PROBABLY MAKE THIS RUN FASTER IN GENERAL,
        % ANYWAY
        mean_mrl_acrossSessions = nan(numSessions, numRegions, numEventTypes, length(phase_f), length(amp_f), length(t));
        mean_mrl_z_acrossSessions = nan(numSessions, numRegions, numEventTypes, length(phase_f), length(amp_f), length(t));
        numSessions_perRegion = zeros(1, numRegions);
        for iSession = 1 : numSessions
            fprintf('analyzing %s, session %d of %d\n', sessionList{iSession}, iSession, numSessions)
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir'); continue; end

            region_meanName = ['gaborPhaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_meanName = fullfile(phaseAmp_sessionDir, region_meanName);

            if ~exist(region_meanName, 'file'); continue; end
            load(region_meanName);
            
            sessionRegionList = region_phaseAmp_metadata.regionList;
            numSessionRegions = length(sessionRegionList);
            
            for iSessionRegion = 1 : numSessionRegions
                regionIdx = find(strcmpi(sessionRegionList{iSessionRegion}, ROI_list));
                if isempty(regionIdx); continue; end
                
                fprintf('%s, %d of %d regions\n', sessionRegionList{iSessionRegion}, iSessionRegion, numSessionRegions)
                
                if size(mean_mrl_acrossSessions, 3) == 1    % indexing issue related to "squeezing" too much if there is only one eventType
                    mean_mrl_acrossSessions(iSession, regionIdx, 1, :, :, :) = ...
                        squeeze(mean_mrl(iSessionRegion, :, :, :, :));
                    mean_mrl_z_acrossSessions(iSession, regionIdx, 1, :, :, :) = ...
                        squeeze(mean_mrl_z(iSessionRegion, :, :, :, :));
                else
                    mean_mrl_acrossSessions(iSession, regionIdx, :, :, :, :) = ...
                        squeeze(mean_mrl(iSessionRegion, :, :, :, :));
                    mean_mrl_z_acrossSessions(iSession, regionIdx, :, :, :, :) = ...
                        squeeze(mean_mrl_z(iSessionRegion, :, :, :, :));
                end
                
                numSessions_perRegion(regionIdx) = numSessions_perRegion(regionIdx) + 1;

            end
            
        end    % for iSession...
       %% 
        mean_mrl_byRegion = squeeze(nanmean(mean_mrl_acrossSessions, 1));
        mean_mrl_z_byRegion = squeeze(nanmean(mean_mrl_z_acrossSessions, 1));
        
        phaseAmpSummary_metadata.numSessions_perRegion = numSessions_perRegion;
        phaseAmpSummary_metadata.numSessionsComplete = iSession;
        save(mat_saveName, 'mean_mrl_byRegion','phaseAmpSummary_metadata');
        save(mat_z_saveName, 'mean_mrl_z_byRegion','phaseAmpSummary_metadata');

        %%
        numPages = 0;
        numRegionPlots = 0;
        h_fig = cell(length(var_to_plot), length(plotTypes)); h_axes = cell(length(var_to_plot), length(plotTypes));

        for iRegion = 1 : numRegions
            numRegionPlots = numRegionPlots + 1;
            
            rowNum = rem(numRegionPlots, regions_per_page);
            if rowNum == 1
                for iVar = 1 : length(var_to_plot)
                    for iPlotType = 1 : length(plotTypes)
                        h_fig{iVar, iPlotType} = zeros(1, length(plotFreqs{iPlotType}));
                        h_axes{iVar, iPlotType} = zeros(length(plotFreqs{iPlotType}), figProps.m, figProps.n);
                        for iFreq = 1 : length(plotFreqs{iPlotType})
                            [h_fig{iVar, iPlotType}(iFreq), h_axes{iVar, iPlotType}(iFreq, :, :)] = createFigPanels5(figProps);
                        end
                    end
                end
                page_regionList = ROI_list{iRegion};
                numSessions_perRegionList = num2str(numSessions_perRegion(iRegion));
%                         page_locList = ch.location.subclass;
                numPages = numPages + 1;
            else
                page_regionList = [page_regionList ', ' ROI_list{iRegion}];
                numSessions_perRegionList = [numSessions_perRegionList ', ' num2str(numSessions_perRegion(iRegion))];
%                         page_locList = [page_locList ', ' ch.location.subclass];
            end
            %%
            if rowNum == 0; rowNum = regions_per_page; end

            for iVar = 1 : length(var_to_plot)
                for iPlotType = 1 : length(plotTypes)
                    for iFreq = 1 : length(plotFreqs{iPlotType})
                        for iEventType = 1 : numEventTypes
                            switch plotTypes{iPlotType}
                                case 'const_phase_f',
                                    if strcmpi(var_to_plot{iVar}, 'mrl')
                                        toPlot = squeeze(mean_mrl_byRegion(iRegion, iEventType, phase_freq_idx(iFreq), :, :));
                                        colorLim = mrl_clim;
                                        textStr{1} = sprintf('mrl, phase-amplitude coupling, constant phase freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                    else
                                        toPlot = squeeze(mean_mrl_z_byRegion(iRegion, iEventType, phase_freq_idx(iFreq), :, :));
                                        colorLim = z_clim;
                                        textStr{1} = sprintf('mrl z-score, phase-amplitude coupling, constant phase freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                    end
                                    x = t;
                                    y = 1:length(amp_f);
                                    y_ticks = amp_freqTick_idx;
                                    x_ticks = t_ticks;
                                    yticklabel = amp_freqTick_label;
                                    xticklabel = x_ticks;
                                case 'const_amp_f',
                                    if strcmpi(var_to_plot{iVar}, 'mrl')
                                        toPlot = squeeze(mean_mrl_byRegion(iRegion, iEventType, :, amp_freq_idx(iFreq), :));
                                        colorLim = mrl_clim;
                                        textStr{1} = sprintf('mrl, phase-amplitude coupling, constant amplitude freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                    else
                                        toPlot = squeeze(mean_mrl_z_byRegion(iRegion, iEventType, :, amp_freq_idx(iFreq), :));
                                        colorLim = z_clim;
                                        textStr{1} = sprintf('mrl z-score, phase-amplitude coupling, constant amplitude freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                    end
                                    x = t;
                                    y = 1:length(phase_f);
                                    y_ticks = phase_freqTick_idx;
                                    x_ticks = t_ticks;
                                    yticklabel = phase_freqTick_label;
                                    xticklabel = x_ticks;
                                case 'averaged_t',
                                    if strcmpi(var_to_plot{iVar}, 'mrl')
                                        toPlot = squeeze(mean(mean_mrl_byRegion(iRegion, iEventType, :, :, :), 5))';
                                        colorLim = mrl_clim;
                                        textStr{1} = 'mrl, phase-amplitude coupling, average across time';
                                    else
                                        toPlot = squeeze(mean(mean_mrl_z_byRegion(iRegion, iEventType, :, :, :), 5))';
                                        colorLim = z_clim;
                                        textStr{1} = 'mrl z-score, phase-amplitude coupling, average across time';
                                    end
                                    x = 1:length(phase_f);
                                    y = 1:length(amp_f);
                                    x_ticks = phase_freqTick_idx;
                                    y_ticks = amp_freqTick_idx;
                                    yticklabel = amp_freqTick_label;
                                    xticklabel = phase_freqTick_label;
                            end

                            axes(h_axes{iVar, iPlotType}(iFreq, rowNum, iEventType));
                            imagesc(x,y,toPlot);    % need to check that toPlot is in the correct orientation
                            set(gca,'ydir','normal',...
                                    'clim',colorLim,...
                                    'xtick',x_ticks,...
                                    'ytick',y_ticks);
                            colormap jet;
                                        
                            if rowNum == 1
                                title(region_phaseAmp_metadata.eventList{iEventType});
                            end
                            if rowNum < regions_per_page
                                set(gca,'xticklabel',[]);
                            else
                                set(gca,'xticklabel',xticklabel);
                            end
                            if iEventType > 1
                                set(gca,'yticklabel',[]);
                            else
                                set(gca,'yticklabel',yticklabel);
                            end

                        end    % for iEventType...
                        
                        if (rowNum == regions_per_page || numRegionPlots == numRegions)
                            h_figAxes = createFigAxes(h_fig{iVar, iPlotType}(iFreq));
                            axes(h_figAxes);

                            textStr{2} = ['Trial type: ' region_phaseAmp_metadata.trialType];
                            textStr{3} = page_regionList;
                            textStr{4} = ['Number of sessions in average: ' numSessions_perRegionList];
                            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                            if numPages == 1
                                export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-nocrop');
                            else
                                export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-append','-nocrop');
                            end
                            close(h_fig{iVar, iPlotType}(iFreq));
                        end
                            
                    end    % for iFreq...
                end    % for iPlotType...

            end    % for iVar...

        end    % for iRegion...
            
    end    % for iTrialType...
            
end    % for i_chDB
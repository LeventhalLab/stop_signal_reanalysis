% script_plot_phaseAmp_Gabors

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed_Gabors';

makePlots = true;
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};
numTrialTypes = length(trialTypeList);

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
channels_per_page = 5;
z_clim = [-4 4];
mrl_clim = [0 1e-2];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

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

cmap = 'jet';

if exist('region_phaseAmp_metadata','var')
    clear region_phaseAmp_metadata;
end
for i_chDB = 1 : 1%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    end
    
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        continue;
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iTrialType = 2 : 2%length(trialTypeList)
        trialType = trialTypeList{iTrialType};
        
%         switch iTrialType
%             case 1,
%                 sessions_to_plot = [];
%             case 2,
%                 sessions_to_plot = [];
%             case 3,
%                 sessions_to_plot = [];
%             case 4,
%                 sessions_to_plot = [];
%             case 5,
%                 sessions_to_plot = [9,11,17];
%             case 6,
%                 sessions_to_plot = [32];
%             case 7,
%                 sessions_to_plot = [];
%         end
        
        for iCh = 1 : length(channels)
            ch = channels{iCh};
            session_paDir = fullfile(subject_phaseAmpdir,sessionList{1},[sessionList{1} '_' trialType]);
            test_paName = fullfile(session_paDir,[ch.name '_' trialType '_phase_amp_Gabor.mat']);
            if ~exist(test_paName,'file');continue;end
            break;
        end
        load(test_paName);
        t = phaseAmp_metadata.t; f = phaseAmp_metadata.f;
        phase_f = phaseAmp_metadata.phase_f; amp_f = phaseAmp_metadata.amp_f;
        numEvents = length(phaseAmp_metadata.eventList);
        numSamps  = length(t);
        
        for iSession = 1 : numSessions
            
%             iSession = sessions_to_plot(iPlotSession);
            sprintf('%d of %d sessions, %s', iSession, numSessions, sessionList{iSession})
            
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir'); continue; end
    
            fig_saveName = cell(length(var_to_plot), length(plotTypes));
            for iVar = 1 : length(var_to_plot)
                for iPlotType = 1 : length(plotTypes)
                    fig_saveName{iVar, iPlotType} = cell(1,length(plotFreqs{iPlotType}));
                        
                    for iFreq = 1 : length(plotFreqs{iPlotType})
                        freqStr = sprintf('%04.0f',plotFreqs{iPlotType}(iFreq));
                        fig_saveName{iVar, iPlotType}{iFreq} = ...
                            ['gaborPhaseAmp_' trialTypeList{iTrialType} '_' var_to_plot{iVar} '_' plotTypes{iPlotType} '_' freqStr '_' sessionList{iSession}];
                        fig_saveName{iVar, iPlotType}{iFreq} = fullfile(phaseAmp_sessionDir, fig_saveName{iVar, iPlotType}{iFreq});
                    end

                end    % for iType...
            end
            if exist([fig_saveName{1,1}{1} '.pdf'], 'file'); continue; end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.locationSubClass = ROI_list;
        
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
        
            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels); continue; end

%             regionList = getRegionsfromChannelDB(sessionChannels);
%             numRegions = length(regionList);
        
            totalSessionChannels = length(sessionChannels);

            region_mean_saveName = ['gaborPhaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_mean_saveName = fullfile(phaseAmp_sessionDir, region_mean_saveName);

            if exist(region_mean_saveName, 'file'); continue; end  % comment this line out to make plots whether means have been saved or not
        
            phaseAmp_surr_name_mrl = [sessionChannels{1}.name '_' trialType '_phase_amp_Gabor_surrogates.mat'];
            phaseAmp_surr_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_surr_name_mrl);
            if ~exist(phaseAmp_surr_name_mrl, 'file'); continue; end
                
            phaseAmp_name_mrl = [sessionChannels{1}.name '_' trialType '_phase_amp_Gabor.mat'];
            phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
            if ~exist(phaseAmp_name_mrl, 'file'); continue; end
            load( phaseAmp_name_mrl );
            
            numEventTypes = length(phaseAmp_metadata.eventList);
            Fs = phaseAmp_metadata.Fs;
            numSamps = size(mrv, 4);

            tmin = phaseAmp_metadata.eventtWin(1);
            eventtWin = phaseAmp_metadata.eventtWin;

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
            t_ticks = [phaseAmp_metadata.eventtWin(1),0,phaseAmp_metadata.eventtWin(2)];
            
            figProps.n = numEventTypes;
            figProps.colSpacing = ones(1, figProps.n) * 0.5;
            figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                         2 * sideMargins - ...
                                                         sum(figProps.colSpacing)) / figProps.n;
            numPages = 0;
            mrl_by_region = cell(1, numRegions);
            mrl_z_by_region = cell(1, numRegions);
            mean_mrl = zeros(numRegions, numEventTypes, length(phase_f), length(amp_f), length(t));
            mean_mrl_z = zeros(numRegions, numEventTypes, length(phase_f), length(amp_f), length(t));
            
            num_ch_per_region = zeros(1, numRegions);
            numChPlots = 0;
            h_fig = cell(length(var_to_plot), length(plotTypes)); h_axes = cell(length(var_to_plot), length(plotTypes));
            for iRegion = 1 : numRegions
            
                cp = initChanParams();
                cp.locationSubClass = ROI_list{iRegion};            
                chList = extractChannels(cp, sessionChannels);
                if isempty(chList); continue; end
                regionChannels = sessionChannels(chList);
                numCh = length(regionChannels);
           
                mrl_by_region{iRegion} = zeros(numCh, numEventTypes, length(phase_f), length(amp_f), length(t));
                mrl_z_by_region{iRegion} = zeros(numCh, numEventTypes, length(phase_f), length(amp_f), length(t));
                num_ch_per_region(iRegion) = numCh;
            
                for iCh = 1 : numCh
    %               iCh
                    ch = regionChannels{iCh};
                    phaseAmp_name_mrl = [ch.name '_' trialType '_phase_amp_Gabor_surrogates.mat'];
                    phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
                    if ~exist(phaseAmp_name_mrl, 'file'); continue; end
                    
                    phaseAmp_surr_name_mrl = [ch.name '_' trialType '_phase_amp_Gabor_surrogates.mat'];
                    phaseAmp_surr_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_surr_name_mrl);
                    if ~exist(phaseAmp_surr_name_mrl, 'file'); continue; end
            
                    load( phaseAmp_name_mrl );
                    load( phaseAmp_surr_name_mrl );
%                     mrv = re_mrv + 1i*im_mrv;
                    mrv_z = (abs(mrv) - surrogate_mean) ./ surrogate_std;

                    numChPlots = numChPlots + 1;
                    rowNum = rem(numChPlots, channels_per_page);
                    if rowNum == 1
                        if makePlots
                            for iVar = 1 : length(var_to_plot)
                                for iPlotType = 1 : length(plotTypes)
                                    h_fig{iVar, iPlotType} = zeros(1, length(plotFreqs{iPlotType}));
                                    h_axes{iVar, iPlotType} = zeros(length(plotFreqs{iPlotType}), figProps.m, figProps.n);
                                    for iFreq = 1 : length(plotFreqs{iPlotType})
                                        [h_fig{iVar, iPlotType}(iFreq), h_axes{iVar, iPlotType}(iFreq, :, :)] = createFigPanels5(figProps);
                                    end
                                end
                            end
                        end
                        page_chList = ch.name;
                        page_locList = ch.location.subclass;
                        numPages = numPages + 1;
                    else
                        page_chList = [page_chList ', ' ch.name];
                        page_locList = [page_locList ', ' ch.location.subclass];
                    end
                    if rowNum == 0; rowNum = channels_per_page; end

                    mrl_by_region{iRegion}(iCh, :, :, :, :) = abs(mrv);
                    mrl_z_by_region{iRegion}(iCh, :, :, :, :) = (abs(mrv) - surrogate_mean) ./ surrogate_std;
                
                    if makePlots                             
                        
                        for iVar = 1 : length(var_to_plot)
                            for iPlotType = 1 : length(plotTypes)
                                for iFreq = 1 : length(plotFreqs{iPlotType})
                                    for iEventType = 1 : numEventTypes
                                        switch plotTypes{iPlotType}
                                            case 'const_phase_f',
                                                if strcmpi(var_to_plot{iVar}, 'mrl')
                                                    toPlot = squeeze(mrl_by_region{iRegion}(iCh, iEventType, phase_freq_idx(iFreq), :, :));
                                                    colorLim = mrl_clim;
                                                    textStr{1} = sprintf('mrl, phase-amplitude coupling, constant phase freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                                else
                                                    toPlot = squeeze(mrl_z_by_region{iRegion}(iCh, iEventType, phase_freq_idx(iFreq), :, :));
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
                                                    toPlot = squeeze(mrl_by_region{iRegion}(iCh, iEventType, :, amp_freq_idx(iFreq), :));
                                                    colorLim = mrl_clim;
                                                    textStr{1} = sprintf('mrl, phase-amplitude coupling, constant amplitude freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                                else
                                                    toPlot = squeeze(mrl_z_by_region{iRegion}(iCh, iEventType, :, amp_freq_idx(iFreq), :));
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
                                                    toPlot = squeeze(mean(mrl_by_region{iRegion}(iCh, iEventType, :, :, :), 5))';
                                                    colorLim = mrl_clim;
                                                    textStr{1} = 'mrl, phase-amplitude coupling, average across time';
                                                else
                                                    toPlot = squeeze(mean(mrl_z_by_region{iRegion}(iCh, iEventType, :, :, :), 5))';
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
                                            title(phaseAmp_metadata.eventList{iEventType});
                                        end
                                        if rowNum < channels_per_page
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

                                    if (rowNum == channels_per_page || numChPlots == totalSessionChannels)
                                        h_figAxes = createFigAxes(h_fig{iVar, iPlotType}(iFreq));
                                        axes(h_figAxes);

                                        textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                                        textStr{3} = page_chList;
                                        textStr{4} = page_locList;
                                        textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                                    
                                        pdfSaveName = [fig_saveName{iVar, iPlotType}{iFreq} '_' num2str(numPages) '.pdf'];
%                                         if numPages == 1
%                                             export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-nocrop');
%                                         else
%                                             export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-append','-nocrop');
%                                         end
                                        export_fig(pdfSaveName,'-pdf','-q101','-painters','-nocrop');
                                        close(h_fig{iVar, iPlotType}(iFreq));
                                    end

                                end    % for iFreq...
                            end    % for iPlotType...
                        end    % for iVar...

                    end    % if makePlots

                end    % for iCh...
                % create averages within individual brain regions for each session
                mean_mrl(iRegion, : ,:, :, :) = mean(mrl_by_region{iRegion}, 1);
                mean_mrl_z(iRegion, :, :, :, :) = mean(mrl_z_by_region{iRegion}, 1);
            
            end    % for iRegion...
            
            region_phaseAmp_metadata = phaseAmp_metadata;
            region_phaseAmp_metadata.regionList = ROI_list;
            region_phaseAmp_metadata.num_ch_per_region = num_ch_per_region;
            region_phaseAmp_metadata.phase_f = phase_f;
            region_phaseAmp_metadata.amp_f = amp_f;
            save(region_mean_saveName, 'mean_mrl', 'mean_mrl_z', 'region_phaseAmp_metadata');
%             save(region_mean_saveName, 'mean_mrl', 'region_phaseAmp_metadata');
            
            % create a new page of figures to show the region means for
            % each session
            if makePlots
                for iRegion = 1 : numRegions
                    
                    rowNum = rem(iRegion, channels_per_page);

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
                        page_numChList  = ['Number of Channels per Region: ' num2str(num_ch_per_region(iRegion))];
                        numPages = numPages + 1;
                    else
                        page_regionList = [page_regionList ', ' ROI_list{iRegion}];
                        page_numChList  = [page_numChList ', ' num2str(num_ch_per_region(iRegion))];
                    end
                    if rowNum == 0; rowNum = channels_per_page; end
                        
                    for iVar = 1 : length(var_to_plot)
                        for iPlotType = 1 : length(plotTypes)
                            for iFreq = 1 : length(plotFreqs{iPlotType})
                                for iEventType = 1 : numEventTypes
                                    switch plotTypes{iPlotType}

                                        case 'const_phase_f',
                                            if strcmpi(var_to_plot{iVar}, 'mrl')
                                                toPlot = squeeze(mean_mrl(iRegion, iEventType, phase_freq_idx(iFreq), :, :));
                                                colorLim = mrl_clim;
                                                textStr{1} = sprintf('mrl, phase-amplitude coupling, constant phase freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                            else
                                                toPlot = squeeze(mean_mrl_z(iRegion, iEventType, phase_freq_idx(iFreq), :, :));
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
                                                toPlot = squeeze(mean_mrl(iRegion, iEventType, :, amp_freq_idx(iFreq), :));
                                                colorLim = mrl_clim;
                                                textStr{1} = sprintf('mrl, phase-amplitude coupling, constant amplitude freq = %f Hz', plotFreqs{iPlotType}(iFreq));
                                            else
                                                toPlot = squeeze(mean_mrl_z(iRegion, iEventType, :, amp_freq_idx(iFreq), :));
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
                                                toPlot = squeeze(mean(mean_mrl(iRegion, iEventType, :, :, :), 5))';
                                                colorLim = mrl_clim;
                                                textStr{1} = 'mrl, phase-amplitude coupling, average across time';
                                            else
                                                toPlot = squeeze(mean(mean_mrl_z(iRegion, iEventType, :, :, :), 5))';
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
                                        title(phaseAmp_metadata.eventList{iEventType});
                                    end
                                    if rowNum < channels_per_page
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

                                if (rowNum == channels_per_page || iRegion == numRegions)
                                    h_figAxes = createFigAxes(h_fig{iVar, iPlotType}(iFreq));
                                    axes(h_figAxes);

                                    textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                                    textStr{3} = page_regionList;
                                    textStr{4} = page_numChList;
                                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                                    pdfSaveName = [fig_saveName{iVar, iPlotType}{iFreq} '_' num2str(numPages) '.pdf'];
                                    export_fig(pdfSaveName,'-pdf','-q101','-painters','-append','-nocrop');
                                    close(h_fig{iVar, iPlotType}(iFreq));
                                end

                            end    % for iFreq...  
                        end    % for iPLotType
                    end    % for iVar...
                end    % for iRegion...
                %%
                % now merge pdfs for each plot type into a single file
                % THERE'S A PROBLEM WITH THIS - THE IMAGES SHOW UP
                % PIXELATED IN THE COMPOSITE FILES, BUT LOOK FINE IN THE
                % ORIGINAL FILES
%                 for iVar = 1 : length(var_to_plot)
%                     for iPlotType = 1 : length(plotTypes)
%                         for iFreq = 1 : length(plotFreqs{iPlotType})
%                             full_pdfSaveName = [fig_saveName{iVar, iPlotType}{iFreq} '.pdf'];
%                             for iPage = 1 : numPages
%                                 pdf_to_append = [fig_saveName{iVar, iPlotType}{iFreq} '_' num2str(iPage) '.pdf'];
%                                 
%                                 append_pdfs(full_pdfSaveName, pdf_to_append);
%                                 delete(pdf_to_append);
%                             end
%                         end
%                     end
%                 end
                                
                            
            end    % if makePlots...
                    
        end    % end for iSession...
        
    end    % for iTrialType...
    
end    % for i_chDB...
                
            
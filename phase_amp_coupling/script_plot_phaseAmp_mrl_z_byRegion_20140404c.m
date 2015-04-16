% script_plot_phaseAmp_mrl_byRegion

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';

makePlots = true;
eventList{1} = {'noseCenterIn'};
eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
eventList{3} = eventList{2};
eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any'};
numTrialTypes = length(trialTypeList);
numPlotTypes  = 3;

if exist('fig_saveName','var');clear('fig_saveName');end

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
channels_per_page = 5;
mrl_colorLim = [0 5e-3];
z_colorLim = [-4 4];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

for i_chDB = 3 : 3%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iTrialType = 1 : 1%length(trialTypeList)
        trialType = trialTypeList{iTrialType};

        for iSession = 1 : numSessions
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession});
            if ~exist(phaseAmp_sessionDir, 'dir'); continue; end
            
            hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
            hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
            
            metadata_filename_1Hz   = [sessionList{iSession} 'hilbert_metadata.mat'];
            metadata_filename_025Hz = [sessionList{iSession} 'hilbert_metadata.mat'];
            metadata_filename_1Hz   = fullfile(hilbert_sessionDir_1Hz, metadata_filename_1Hz);
            metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename_025Hz);
            
            if ~exist(metadata_filename_1Hz, 'file'); continue; end
            
            md_1Hz   = load(metadata_filename_1Hz);
            md_025Hz = load(metadata_filename_025Hz);
            
            centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
            centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
            centerFreqs = [centerFreqs_025Hz; centerFreqs_1Hz];

            fig_saveName{1} = ['mrl_phaseAmpPlots_' trialTypeList{iTrialType} '_' sessionList{iSession} '.pdf'];
            fig_saveName{1} = fullfile(phaseAmp_sessionDir, fig_saveName{1});
%             if exist(fig_saveName{1}, 'file'); continue; end
            
            fig_saveName{2} = ['mrl_z_phaseAmpPlots_' trialTypeList{iTrialType} '_' sessionList{iSession} '.pdf'];
            fig_saveName{2} = fullfile(phaseAmp_sessionDir, fig_saveName{2});
            
            fig_saveName{3} = ['mrl_trialz_phaseAmpPlots_' trialTypeList{iTrialType} '_' sessionList{iSession} '.pdf'];
            fig_saveName{3} = fullfile(phaseAmp_sessionDir, fig_saveName{3});
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
        
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);
        
            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
            sessionChannels = excludeChannels(cp, sessionChannels);

            regionList = getRegionsfromChannelDB(sessionChannels);
            numRegions = length(regionList);
        
            totalSessionChannels = length(sessionChannels);

            region_mean_saveName = ['mrl_phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_mean_saveName = fullfile(phaseAmp_sessionDir, region_mean_saveName);
            region_z_mean_saveName = ['mrl_z_phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_z_mean_saveName = fullfile(phaseAmp_sessionDir, region_z_mean_saveName);
            region_trialz_mean_saveName = ['mrl_trialz_phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_trialz_mean_saveName = fullfile(phaseAmp_sessionDir, region_trialz_mean_saveName);

            if exist(region_mean_saveName, 'file'); continue; end  % comment this line out to make plots whether means have been saved or not
        
            phaseAmp_name = [sessionChannels{1}.name '_' trialType '_phase_amp_mrl.mat'];
            surr_phaseAmp_name = [sessionChannels{1}.name '_surrTrial_t_' trialType '_phase_amp_mrl.mat'];
            trial_surr_phaseAmp_name = [sessionChannels{1}.name '_surrogate_phase_amp_mrl.mat'];
            phaseAmp_name = fullfile(phaseAmp_sessionDir, phaseAmp_name);
            surr_phaseAmp_name = fullfile(phaseAmp_sessionDir, surr_phaseAmp_name);
            trial_surr_phaseAmp_name = fullfile(phaseAmp_sessionDir, trial_surr_phaseAmp_name);
            if ~exist(phaseAmp_name, 'file'); continue; end
            if ~exist(surr_phaseAmp_name, 'file'); continue; end
            if ~exist(trial_surr_phaseAmp_name, 'file'); continue; end
            load( phaseAmp_name );
        
            numEventTypes = length(phaseAmp_metadata.eventList);
            Fs = phaseAmp_metadata.Fs;
            numSamps = size(re_mrv, 4);
            eventtWin = phaseAmp_metadata.eventtWin;
            if length(eventtWin) == 1
                phaseAmp_metadata.eventtWin(2) = eventtWin(1) + phaseAmp_metadata.analysisWin;
                eventtWin = phaseAmp_metadata.eventtWin;
                save( phaseAmp_name, 're_mrv', 'im_mrv', 'phaseAmp_metadata');   % corrects an oversight in the code that wrote out these structures
            end
            tmin = eventtWin(1);
%             analysisWin = phaseAmp_metadata.analysisWin;
%             stepSize = phaseAmp_metadata.stepSize;
            
%             winStartTimes = tmin : stepSize : (eventtWin(2) - analysisWin);
%             numSteps = length(winStartTimes);
            
            low_freq_range = phaseAmp_metadata.low_freq_range;
            high_freq_range = phaseAmp_metadata.high_freq_range;
        
            low_freq_idx  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
            high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

            low_freqs  = centerFreqs_025Hz(low_freq_idx);
            high_freqs = centerFreqs_1Hz(high_freq_idx);
            
            figProps.n = numEventTypes;
            figProps.colSpacing = ones(1, figProps.n) * 0.5;
            figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                         2 * sideMargins - ...
                                                         sum(figProps.colSpacing)) / figProps.n;
            numPages = 0;
            mrl_by_region = cell(1, numRegions);
            mrl_z_by_region = cell(1, numRegions);
            mrl_trialz_by_region = cell(1, numRegions);
            mean_mrl = zeros(numRegions, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
            mean_mrl_z = zeros(numRegions, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
            mean_mrl_trialz = zeros(numRegions, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
            
            num_ch_per_region = zeros(1, numRegions);
            numChPlots = 0;
            for iRegion = 1 : numRegions
            
                cp = initChanParams();
                cp.locationName = regionList{iRegion};            
                chList = extractChannels(cp, sessionChannels);
                if isempty(chList); continue; end
                regionChannels = sessionChannels(chList);
                numCh = length(regionChannels);
           
                mrl_by_region{iRegion} = zeros(numCh, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
                mrl_z_by_region{iRegion} = zeros(numCh, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
                mrl_trialz_by_region{iRegion} = zeros(numCh, numEventTypes, length(low_freqs), length(high_freqs), numSamps);
                mrl_z = zeros(numEventTypes, length(low_freqs), length(high_freqs), numSamps);
                mrl_trialz = zeros(numEventTypes, length(low_freqs), length(high_freqs), numSamps);
                num_ch_per_region(iRegion) = numCh;
            
                for iCh = 1 : numCh
    %               iCh
                    ch = regionChannels{iCh};
                    phaseAmp_name = [ch.name '_' trialType '_phase_amp_mrl.mat'];
                    surr_phaseAmp_name = [ch.name '_surrTrial_t_' trialType '_phase_amp_mrl.mat'];
                    phaseAmp_name = fullfile(phaseAmp_sessionDir, phaseAmp_name);
                    surr_phaseAmp_name = fullfile(phaseAmp_sessionDir, surr_phaseAmp_name);
                    if ~exist(phaseAmp_name, 'file'); continue; end
                    if ~exist(surr_phaseAmp_name, 'file'); continue; end
                    
                    load( phaseAmp_name );
                    mrv = re_mrv + 1i*im_mrv;
                    mrl = abs(mrv);
%                     fullsigSurr = load( surr_phaseAmp_name );
                    trialSurr   = load( surr_phaseAmp_name );
                    
                    mrl_trialz = (mrl - trialSurr.surrogate_mean) ./ trialSurr.surrogate_std;
                    
%                     for iEventType = 1 : numEventTypes
%                         for i_f1 = 1 : length(low_freqs)
%                             for i_f2 = 1 : length(high_freqs)
%                                 for i_t = 1 : numSamps
% %                                     mrl_z(iEventType, i_f1, i_f2, i_t) = ...
% %                                         (abs(mrv(iEventType, i_f1, i_f2, i_t)) - fullsigSurr.surrogate_mean(i_f1, i_f2)) / fullsigSurr.surrogate_std(i_f1, i_f2);
%                                     mrl_trialz(iEventType, i_f1, i_f2, i_t) = ...
%                                         (abs(mrv(iEventType, i_f1, i_f2, i_t)) - trialSurr.surrogate_mean(i_f1, i_f2)) / trialSurr.surrogate_std(i_f1, i_f2);
%                                 end
%                             end
%                         end
%                     end
                    




%%%%%% WORKING HERE...



                    numChPlots = numChPlots + 1;
                    rowNum = rem(numChPlots, channels_per_page);
                    if rowNum == 1
                        if makePlots
                            h_fig = zeros(1,3); h_axes = zeros(3, figProps.m, figProps.n);
                            for iPlotType = 1 : 3
                                [h_fig(iPlotType), h_axes(iPlotType,:)] = createFigPanels5(figProps);
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

                    mrl_by_region{iRegion}(iCh, :, :, :, :) = mrl;
                    mrl_z_by_region{iRegion}(iCh, :, :, :, :) = mrl_z;
                    mrl_trialz_by_region{iRegion}(iCh, :, :, :, :) = mrl_trialz;
                
                    for iEventType = 1 : numEventTypes

%                         toPlot = log10(squeeze(mrl(iEventType, 1, :, :))');
%                         toPlot_z = squeeze(mrl_z(iEventType, 1, :, :))';

                        toPlot = squeeze(mean(mrl(iEventType, :, :, :), 4))';   % average over time
                        toPlot_z = squeeze(mean(mrl_z(iEventType, :, :, :), 4))';
                        toPlot_trialz = squeeze(mean(mrl_trialz(iEventType, :, :, :), 4))';
                        
                        if makePlots
                            for iPlotType = 1 : 3
                                axes(h_axes(iPlotType, rowNum, iEventType));
                                switch iPlotType
                                    case 1,
                                        imagesc(low_freqs, high_freqs, toPlot);
                                        set(gca,'ydir','normal');%,'clim',mrl_colorLim);
                                    case 2,
                                        imagesc(low_freqs, high_freqs, toPlot_z);
                                        set(gca,'ydir','normal');%,'clim',z_colorLim);
                                    case 3,
                                        imagesc(low_freqs, high_freqs, toPlot_trialz);
                                        set(gca,'ydir','normal');%,'clim',z_colorLim);
                                end
                                colorbar
                                
                                if rowNum == 1
                                    title(phaseAmp_metadata.eventList{iEventType});
                                end
                                if rowNum < channels_per_page
                                    set(gca,'xticklabel',[]);
                                end
                                if iEventType > 1
                                    set(gca,'yticklabel',[]);
                                end
                            end    % for iPlotType
                        end
                    end

                    if (rowNum == channels_per_page || numChPlots == totalSessionChannels) && makePlots
                        for iPlotType = 1 : 3
                            h_figAxes = createFigAxes(h_fig(iPlotType));
                            axes(h_figAxes);
                            
                            switch iPlotType
                                case 1,
                                    textStr{1} = sprintf('phase-amplitude coupling, twin = [%d, %d] around %s', eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                                case 2,
                                    textStr{1} = sprintf('phase-amplitude coupling full signal z-score, twin = [%d, %d] around %s', eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                                case 3,
                                    textStr{1} = sprintf('phase-amplitude coupling trials z-score, twin = [%d, %d] around %s', eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                            end
                            textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                            textStr{3} = page_chList;
                            textStr{4} = page_locList;
                            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                            figure(h_fig(iPlotType));
                            
                            if numPages == 1
                                export_fig(fig_saveName{iPlotType},'-pdf','-q101','-painters');
                            else
                                export_fig(fig_saveName{iPlotType},'-pdf','-q101','-painters','-append');
                            end

                            close(h_fig(iPlotType));
                        end    % for iPlotType = 1 : 3
                    end

                end    % for iCh...
                % create averages within individual brain regions for each session
                mean_mrl(iRegion, : ,:, :, :) = mean(mrl_by_region{iRegion}, 1);
                mean_mrl_z(iRegion, : ,:, :, :) = mean(mrl_z_by_region{iRegion}, 1);
                mean_mrl_trialz(iRegion, : ,:, :, :) = mean(mrl_trialz_by_region{iRegion}, 1);
            
            end    % for iRegion...
            region_phaseAmp_metadata = phaseAmp_metadata;
            region_phaseAmp_metadata.regionList = regionList;
            region_phaseAmp_metadata.num_ch_per_region = num_ch_per_region;
            region_phaseAmp_metadata.low_freqs = low_freqs;
            region_phaseAmp_metadata.high_freqs = high_freqs;
            save(region_mean_saveName, 'mean_mrl', 'region_phaseAmp_metadata');
            save(region_z_mean_saveName, 'mean_mrl_z', 'region_phaseAmp_metadata');
        
            % create a new page of figures to show the region means for
            % each session
            if makePlots
                for iRegion = 1 : numRegions
                    rowNum = rem(iRegion, channels_per_page);
                    if rowNum == 1
                        for iPlotType = 1 : 2
                            [h_fig(iPlotType), h_axes(iPlotType,:,:)] = createFigPanels5(figProps);
                        end
                        page_regionList = regionList{iRegion};
                        page_numChList  = ['Number of Channels per Region: ' num2str(num_ch_per_region(iRegion))];
                    else
                        page_regionList = [page_regionList ', ' regionList{iRegion}];
                        page_numChList  = [page_numChList ', ' num2str(num_ch_per_region(iRegion))];
                    end
                    if rowNum == 0; rowNum = channels_per_page; end

                    for iEventType = 1 : numEventTypes
                        toPlot = squeeze(mean(mean_mrl(iRegion, iEventType, :, :, :), 5))';
                        toPlot_z = squeeze(mean(mean_mrl_z(iRegion, iEventType, :, :, :), 5))';
                        toPlot_trialz = squeeze(mean(mean_mrl_trialz(iRegion, iEventType, :, :, :), 5))';

                        for iPlotType = 1 : numPlotTypes
                            axes(h_axes(iPlotType, rowNum, iEventType));
                            switch iPlotType
                                case 1,
                                    imagesc(low_freqs, high_freqs, toPlot);
                                    set(gca,'ydir','normal');%,'clim',mrl_colorLim);
                                case 2,
                                    imagesc(low_freqs, high_freqs, toPlot_z);
                                    set(gca,'ydir','normal');%,'clim',z_colorLim);
                                case 3,
                                    imagesc(low_freqs, high_freqs, toPlot_trialz);
                                    set(gca,'ydir','normal');%,'clim',z_colorLim);
                            end
                            colorbar
                            if rowNum == 1
                                title(phaseAmp_metadata.eventList{iEventType});
                            end
                            if rowNum < channels_per_page
                                set(gca,'xticklabel',[]);
                            end
                            if iEventType > 1
                                set(gca,'yticklabel',[]);
                            end
                        end
                    end

                    if (rowNum == channels_per_page || iRegion == numRegions)
                        for iPlotType = 1 : numPlotType
                            h_figAxes = createFigAxes(h_fig);
                            figure(h_fig(iPlotType));
                            axes(h_figAxes);

                            switch iPlotType
                                case 1
                                    textStr{1} = sprintf('%s phase-amplitude coupling within regions, twin = [%d, %d] around %s', implantID, eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                                case 2,
                                    textStr{1} = sprintf('%s phase-amplitude full signal z-score within regions, twin = [%d, %d] around %s', implantID, eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                                case 3,
                                    textStr{1} = sprintf('%s phase-amplitude trials z-score within regions, twin = [%d, %d] around %s', implantID, eventtWin(1), eventtWin(2), phaseAmp_metadata.eventList{1});
                            end
                            textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                            textStr{3} = page_regionList;
                            textStr{4} = page_numChList;
                            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                            export_fig(fig_saveName{iPlotType},'-pdf','-q101','-painters','-append');

                            close(h_fig(iPlotType));
                        end
                    end
                end    % for iRegion...
            end    
        end    % for iSession...
                    
    end    % end 
    
end    % for i_chDB...
                
            
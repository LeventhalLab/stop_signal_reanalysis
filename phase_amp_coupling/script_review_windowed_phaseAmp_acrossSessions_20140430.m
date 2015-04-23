% script_review_windowed_phaseAmp_acrossSessions

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
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};
numTrialTypes = length(trialTypeList);

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
regions_per_page = 5;
z_clim = [-3 3];
mrl_clim = [0 1e-2];
% colorLim = [-4 -1];

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
                                                  
for i_chDB = 4 : 4%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        continue;
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    % exclude EMG, reference channels
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
    channels = excludeChannels(cp, channels);
    
    cp = initChanParams();
    cp.locationSubClass = {'e2', 'e3', 'e02','e03'};
    channels = excludeChannels(cp, channels);
    if isempty(channels);continue;end
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    allRegionList = getRegionsfromChannelDB( channels );
    num_allRegions = length(allRegionList);
    
    for iTrialType = 5 : 5%length(trialTypeList)
        trialType = trialTypeList{iTrialType};

        fig_saveName = cell(length(var_to_plot), length(plotTypes));
        mat_saveName = [implantID '_phaseAmp_' trialTypeList{iTrialType} '_summary.mat'];
        mat_saveName = fullfile(subject_phaseAmpdir, mat_saveName);
        
        mat_z_saveName = [implantID '_phaseAmp_' trialTypeList{iTrialType} '_z_summary.mat'];
        mat_z_saveName = fullfile(subject_phaseAmpdir, mat_z_saveName);
        
        startSession = 1;
        if exist(mat_saveName,'file')
            load(mat_saveName);
            if phaseAmpSummary_metadata.numSessionsComplete < numSessions    % all the calculations have not been done yet
                startSession = iSession;
            
            else    % all the calculations are already done
                continue;
            end
        
        end
        
        % load the first session mean mrl file to figure out how to
        % allocate memory
        for iSession = 1 : numSessions
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            region_meanName = ['phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_meanName = fullfile(phaseAmp_sessionDir, region_meanName);
            if ~exist(region_meanName, 'file'); continue; end
            load(region_meanName);
            break;
        end
        
        
        for iVar = 1 : length(var_to_plot)
            for iPlotType = 1 : length(plotTypes)
                fig_saveName{iVar, iPlotType} = cell(1,length(plotFreqs{iPlotType}));

                for iFreq = 1 : length(plotFreqs{iPlotType})
                    freqStr = sprintf('%04.0f',plotFreqs{iPlotType}(iFreq));
                    fig_saveName{iVar, iPlotType}{iFreq} = ...
                        ['phaseAmp_' trialTypeList{iTrialType} '_' var_to_plot{iVar} '_' plotTypes{iPlotType} '_' freqStr '.pdf'];
                    fig_saveName{iVar, iPlotType}{iFreq} = fullfile(subject_phaseAmpdir, fig_saveName{iVar, iPlotType}{iFreq});
                end

            end    % for iPlotType...
        end
        
        phaseAmpSummary_metadata.low_freq_range = region_phaseAmp_metadata.low_freq_range;
        phaseAmpSummary_metadata.high_freq_range = region_phaseAmp_metadata.high_freq_range;
        phaseAmpSummary_metadata.eventList = region_phaseAmp_metadata.eventList;
        phaseAmpSummary_metadata.eventtWin = region_phaseAmp_metadata.eventtWin;
        phaseAmpSummary_metadata.Fs = region_phaseAmp_metadata.Fs;
        phaseAmpSummary_metadata.low_freqs = region_phaseAmp_metadata.low_freqs;
        phaseAmpSummary_metadata.high_freqs = region_phaseAmp_metadata.high_freqs;
        phaseAmpSummary_metadata.regionList = allRegionList;
        
        eventtWin = region_phaseAmp_metadata.eventtWin;
        numSamps = size(mean_mrl, 5);
        t = linspace(eventtWin(1), eventtWin(2), numSamps);
        
        low_freqs = region_phaseAmp_metadata.low_freqs;
        high_freqs = region_phaseAmp_metadata.high_freqs;
        
        numEventTypes = length(region_phaseAmp_metadata.eventList);
        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
                                                 
        phase_freq_idx = zeros(size(phase_freq));
        amp_freq_idx = zeros(size(amp_freq));
        for iFreq = 1 : length(phase_freq)
            phase_freq_idx(iFreq) = find(low_freqs == phase_freq(iFreq));
        end
        for iFreq = 1 : length(amp_freq)
            amp_freq_idx(iFreq) = find(high_freqs == amp_freq(iFreq));
        end
        
%         mean_mrl_acrossSessions = nan(numSessions, num_allRegions, numEventTypes, length(low_freqs), length(high_freqs), length(t));
%         mean_mrl_z_acrossSessions = nan(numSessions, num_allRegions, numEventTypes, length(low_freqs), length(high_freqs), length(t));

        if startSession == 1
            mean_mrl_byRegion = zeros(num_allRegions, numEventTypes, length(low_freqs), length(high_freqs), length(t));
            mean_mrl_z_byRegion = zeros(num_allRegions, numEventTypes, length(low_freqs), length(high_freqs), length(t));

            numSessions_perRegion = zeros(1, num_allRegions);
        end
        for iSession = startSession : numSessions
            disp(sprintf('analyzing %s, session %d of %d', sessionList{iSession}, iSession, numSessions))
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir'); continue; end
            
%             hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
%             hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
            
%             cp = initChanParams();
%             cp.session = sessionList{iSession};
%         
%             session_chList = extractChannels( cp, channels );
%             sessionChannels = channels(session_chList);
        
%             % exclude EMG, reference channels
%             cp = initChanParams();
%             cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
%             sessionChannels = excludeChannels(cp, sessionChannels);

            region_meanName = ['phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_meanName = fullfile(phaseAmp_sessionDir, region_meanName);

            if ~exist(region_meanName, 'file'); continue; end  % comment this line out to make plots whether means have been saved or not
            load(region_meanName);
            
            sessionRegionList = region_phaseAmp_metadata.regionList;
            numSessionRegions = length(sessionRegionList);
            
            for iSessionRegion = 1 : numSessionRegions
                regionIdx = find(strcmpi(sessionRegionList{iSessionRegion}, allRegionList));
                if isempty(regionIdx); continue; end
                
                disp(sprintf('%s, %d of %d regions', sessionRegionList{iSessionRegion}, iSessionRegion, numSessionRegions))
                
%                 if size(mean_mrl_acrossSessions, 3) == 1    % indexing issue related to "squeezing" too much if there is only one eventType
                    mean_mrl_byRegion(regionIdx, :, :, :, :) = ...
                        mean_mrl_byRegion(regionIdx, :, :, :, :) + ...
                        mean_mrl(iSessionRegion, :, :, :, :);
                    mean_mrl_z_byRegion(regionIdx, :, :, :, :) = ...
                        mean_mrl_z_byRegion(regionIdx, :, :, :, :) + ...
                        mean_mrl_z(iSessionRegion, :, :, :, :);
%                 else
%                     mean_mrl_byRegion(iSession, regionIdx, :, :, :, :) = ...
%                         squeeze(mean_mrl(iSessionRegion, :, :, :, :));
%                     mean_mrl_z_byRegion(iSession, regionIdx, :, :, :, :) = ...
%                         squeeze(mean_mrl_z(iSessionRegion, :, :, :, :));
%                 end
                
                numSessions_perRegion(regionIdx) = numSessions_perRegion(regionIdx) + 1;

            end
            
            phaseAmpSummary_metadata.numSessions_perRegion = numSessions_perRegion;
            phaseAmpSummary_metadata.numSessionsComplete = iSession;
            save(mat_saveName, 'mean_mrl_byRegion','phaseAmpSummary_metadata');
            save(mat_z_saveName, 'mean_mrl_z_byRegion','phaseAmpSummary_metadata');
            
            
        end    % for iSession...
        
        for iRegion = 1 : num_allRegions
            mean_mrl_byRegion(iRegion, :, :, :, :) = mean_mrl_byRegion(iRegion, :, :, :, :) / numSessions_perRegion(iRegion);%squeeze(nanmean(mean_mrl_acrossSessions, 1));
            mean_mrl_z_byRegion(iRegion, :, :, :, :) = mean_mrl_z_byRegion(iRegion, :, :, :, :) / numSessions_perRegion(iRegion);%squeeze(nanmean(mean_mrl_z_acrossSessions, 1));
        end
        
  %%   section for making the plots
        numPages = 0;
        numRegionPlots = 0;
        h_fig = cell(length(var_to_plot), length(plotTypes)); h_axes = cell(length(var_to_plot), length(plotTypes));
        for iRegion = 1 : num_allRegions
            
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
                page_regionList = allRegionList{iRegion};
                numSessions_perRegionList = num2str(numSessions_perRegion(iRegion));
%                         page_locList = ch.location.subclass;
                numPages = numPages + 1;
            else
                page_regionList = [page_regionList ', ' allRegionList{iRegion}];
                numSessions_perRegionList = [numSessions_perRegionList ', ' num2str(numSessions_perRegion(iRegion))];
%                         page_locList = [page_locList ', ' ch.location.subclass];
            end
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
                                        y = high_freqs;
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
                                        y = low_freqs;
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
                                        x = low_freqs;
                                        y = high_freqs;
                                end

                                axes(h_axes{iVar, iPlotType}(iFreq, rowNum, iEventType));
                                imagesc(x,y,toPlot);    % need to check that toPlot is in the correct orientation
                                set(gca,'ydir','normal','clim',colorLim);

                                if rowNum == 1
                                    title(region_phaseAmp_metadata.eventList{iEventType});
                                end
                                if rowNum < regions_per_page
                                    set(gca,'xticklabel',[]);
                                end
                                if iEventType > 1
                                    set(gca,'yticklabel',[]);
                                end
                                
                            end    % for iEventType...

                            if (rowNum == regions_per_page || numRegionPlots == num_allRegions)
                                colorbar
                                h_figAxes = createFigAxes(h_fig{iVar, iPlotType}(iFreq));
                                axes(h_figAxes);

                                textStr{2} = ['Trial type: ' region_phaseAmp_metadata.trialType];
                                textStr{3} = page_regionList;
                                textStr{4} = ['Number of sessions in average: ' numSessions_perRegionList];
                                textStr{5} = ['Color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
                                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                                if numPages == 1
                                    export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters');
                                else
                                    export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-append');
                                end
                                close(h_fig{iVar, iPlotType}(iFreq));
                            end
                                

                        end    % for iFreq...
                    end    % for iPlotType...
                end    % for iVar...

            

        end    % for iRegion...
            

            
    end    % for iTrialType...
            
end    % for i_chDB
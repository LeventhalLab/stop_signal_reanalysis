% script_plot_windowed_phaseAmp_coupling_byRegion

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

phase_freq = 1.5;       % Hz, for the time-frequency plots, can only look at one "phase" frequency at a time

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
channels_per_page = 5;
colorLim = [-4 -1];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
                                                  
for i_chDB = 1 : 1%length(chDB_list)
    
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
    
    for iTrialType = 2 : length(trialTypeList)
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
            
            md_1Hz   = load(metadata_filename_1Hz);
            md_025Hz = load(metadata_filename_025Hz);
            
            centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
            centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
            centerFreqs = [centerFreqs_025Hz; centerFreqs_1Hz];

            fig_saveName = ['phaseAmpPlots_' trialTypeList{iTrialType} '_' sessionList{iSession} '.pdf'];
            fig_saveName = fullfile(phaseAmp_sessionDir, fig_saveName);
%             if exist(fig_saveName, 'file'); continue; end
            
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

            region_mean_saveName = ['phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            region_mean_saveName = fullfile(phaseAmp_sessionDir, region_mean_saveName);

            if exist(region_mean_saveName, 'file'); continue; end  % comment this line out to make plots whether means have been saved or not
        
            phaseAmp_name = [sessionChannels{1}.name '_' trialType '_phase_amp.mat'];
            phaseAmp_name = fullfile(phaseAmp_sessionDir, phaseAmp_name);
            if ~exist(phaseAmp_name, 'file'); continue; end
            load( phaseAmp_name );
        
            numEventTypes = length(phaseAmp_metadata.eventList);
            Fs = phaseAmp_metadata.Fs;
            numSteps = size(MI, 2);
            tmin = phaseAmp_metadata.eventtWin(1);
            eventtWin = phaseAmp_metadata.eventtWin;
            analysisWin = phaseAmp_metadata.analysisWin;
            stepSize = phaseAmp_metadata.stepSize;
            
            low_freq_range = phaseAmp_metadata.low_freq_range;
            high_freq_range = phaseAmp_metadata.high_freq_range;
        
            low_freq_idx  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
            high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

            low_freqs  = centerFreqs_025Hz(low_freq_idx);
            high_freqs = centerFreqs_1Hz(high_freq_idx);
            
            t = linspace(eventtWin(1) + analysisWin/2, eventtWin(2) - analysisWin/2, numSteps);
            
            phase_freq_idx = find(low_freqs == phase_freq);
            
            figProps.n = numEventTypes;
            figProps.colSpacing = ones(1, figProps.n) * 0.5;
            figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                         2 * sideMargins - ...
                                                         sum(figProps.colSpacing)) / figProps.n;
            numPages = 0;
            MI_by_region = cell(1, numRegions);
            mean_MI = zeros(numRegions, numEventTypes, length(t), length(low_freqs), length(high_freqs));
            
            num_ch_per_region = zeros(1, numRegions);
            numChPlots = 0;
            for iRegion = 1 : numRegions
            
                cp = initChanParams();
                cp.locationName = regionList{iRegion};            
                chList = extractChannels(cp, sessionChannels);
                if isempty(chList); continue; end
                regionChannels = sessionChannels(chList);
                numCh = length(regionChannels);
           
                MI_by_region{iRegion} = zeros(numCh, numEventTypes, length(t), length(low_freqs), length(high_freqs));
                num_ch_per_region(iRegion) = numCh;
            
                for iCh = 1 : numCh
    %               iCh
                    ch = regionChannels{iCh};
                    phaseAmp_name = [ch.name '_' trialType '_phase_amp.mat'];
                    phaseAmp_name = fullfile(phaseAmp_sessionDir, phaseAmp_name);
                    if ~exist(phaseAmp_name, 'file'); continue; end
                    
                    load( phaseAmp_name );

                    numChPlots = numChPlots + 1;
                    rowNum = rem(numChPlots, channels_per_page);
                    if rowNum == 1
                        if makePlots
                            [h_fig, h_axes] = createFigPanels5(figProps);
                        end
                        page_chList = ch.name;
                        page_locList = ch.location.subclass;
                        numPages = numPages + 1;
                    else
                        page_chList = [page_chList ', ' ch.name];
                        page_locList = [page_locList ', ' ch.location.subclass];
                    end
                    if rowNum == 0; rowNum = channels_per_page; end

                    MI_by_region{iRegion}(iCh, :, :, :, :) = MI;
                
                    for iEventType = 1 : numEventTypes

                        toPlot = log10(squeeze(MI(iEventType, :, phase_freq_idx, :))');

                        if makePlots
                            axes(h_axes(rowNum, iEventType));
                            imagesc(t, high_freqs, toPlot);
                            set(gca,'ydir','normal','clim',colorLim);
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
                    end    % for iEventType...

                    if (rowNum == channels_per_page || numChPlots == totalSessionChannels) && makePlots
                        h_figAxes = createFigAxes(h_fig);
                        axes(h_figAxes);

                        textStr{1} = ['phase-amplitude coupling, phase frequency = ' num2str(phase_freq)];
                        textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                        textStr{3} = page_chList;
                        textStr{4} = page_locList;
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                        if numPages == 1
                            export_fig(fig_saveName,'-pdf','-q101','-painters');
                        else
                            export_fig(fig_saveName,'-pdf','-q101','-painters','-append');
                        end
                        close(h_fig);
                    end

                end    % for iCh...
                % create averages within individual brain regions for each session
                mean_MI(iRegion, : ,:, :, :) = mean(MI_by_region{iRegion}, 1);
            
            end    % for iRegion...
            region_phaseAmp_metadata = phaseAmp_metadata;
            region_phaseAmp_metadata.regionList = regionList;
            region_phaseAmp_metadata.num_ch_per_region = num_ch_per_region;
            region_phaseAmp_metadata.low_freqs = low_freqs;
            region_phaseAmp_metadata.high_freqs = high_freqs;
            save(region_mean_saveName, 'mean_MI', 'region_phaseAmp_metadata');
        
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
                    if rowNum == 0; rowNum = channels_per_page; end
                
                    for iEventType = 1 : numEventTypes
                        toPlot = squeeze(log10(mean_MI(iRegion, iEventType, :, phase_freq_idx, :)))';
                        axes(h_axes(rowNum, iEventType));
                        imagesc(t, high_freqs, toPlot);
                        set(gca,'ydir','normal','clim',colorLim);
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

                    if (rowNum == channels_per_page || iRegion == numRegions)
                        h_figAxes = createFigAxes(h_fig);
                        axes(h_figAxes);

                        textStr{1} = ['Phase-Amplitude Coupling within Regions, phase frequency = ' num2str(phase_freq)];
                        textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                        textStr{3} = page_regionList;
                        textStr{4} = page_numChList;
                        text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

                        export_fig(fig_saveName,'-pdf','-q101','-painters','-append');
                        close(h_fig);
                    end
                end
            end    % for iRegion...
        end    % if makePlots...
                    
    end    % end for iSession...
    
end    % for i_chDB...
                
            
% script_plot_mean_phaseAmp_coupling_byRegion

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

phase_freq = 1.5;       % Hz, for the time-frequency plots, can only look at one "phase" frequency at a time

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;
regions_per_page = 5;
colorLim = [-7 -4];

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = regions_per_page;    % number of rows

sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

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
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    allRegionList = getRegionsfromChannelDB( channels );
    numSessions = length(sessionList);
    num_allRegions = length(allRegionList);
    
    for iTrialType = 1 : length(trialTypeList)
        trialType = trialTypeList{iTrialType};
        
        fig_saveName = ['mean_phaseAmpPlots_' trialTypeList{iTrialType} '.pdf'];
        fig_saveName = fullfile(subject_phaseAmpdir, fig_saveName);
        
        region_mean_saveName = ['mean_phaseAmp_byRegion_' trialTypeList{iTrialType} '.mat'];
        region_mean_saveName = fullfile(subject_phaseAmpdir, region_mean_saveName);
        if ~exist(region_mean_saveName, 'file'); continue; end
        
        % figure out how big to make the mean_regionMI array
        for iSession = 1 : numSessions
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession});
            
            session_region_mean_name = ['phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            session_region_mean_name = fullfile(phaseAmp_sessionDir, session_region_mean_name);
            if ~exist(session_region_mean_name, 'file'); continue; end
            
            load(session_region_mean_name);
            numEventTypes = length(region_phaseAmp_metadata.eventList);
            numSteps = size(mean_MI, 3);
            num_low_freqs = size(mean_MI, 4);
            num_high_freqs = size(mean_MI, 5);
            
            mean_regionMI = zeros(num_allRegions, numEventTypes, numSteps, num_low_freqs, num_high_freqs);
            
            break;
        end
        num_regionSessions = zeros(1, num_allRegions);
        num_channels_per_regionSession = zeros(numSessions, num_allRegions);
        all_region_phaseAmp_metadata = cell(1, numSessions);
        
        for iSession = 1 : numSessions
%             phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession});
%             if ~exist(phaseAmp_sessionDir, 'dir'); continue; end
            
            session_region_mean_name = ['phaseAmp_byRegion_' trialTypeList{iTrialType} '_' sessionList{iSession} '.mat'];
            session_region_mean_name = fullfile(phaseAmp_sessionDir, session_region_mean_name);
            if exist(session_region_mean_name, 'file'); continue; end   % comment this line out to make plots whether means have been saved or not
            
            load(session_region_mean_name);
            if ~isfield(region_phaseAmp_metadata.low_freqs)   % to fix an oversight in early versions of these files
            
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
                
                low_freq_range = region_phaseAmp_metadata.low_freq_range;
                high_freq_range = region_phaseAmp_metadata.high_freq_range;

                low_freq_idx  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
                high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

                low_freqs  = centerFreqs_025Hz(low_freq_idx);
                high_freqs = centerFreqs_1Hz(high_freq_idx);
                
                region_phaseAmp_metadata.low_freqs = low_freqs;
                region_phaseAmp_metadata.high_freqs = high_freqs;
                
                save(region_mean_saveName, 'mean_MI', 'region_phaseAmp_metadata');
            end
            all_region_phaseAmp_metadata{iSession} = region_phaseAmp_metadata;
            phase_freq_idx = find(low_freqs == phase_freq);

%             cp = initChanParams();
%             cp.session = sessionList{iSession};
%         
%             session_chList = extractChannels( cp, channels );
%             sessionChannels = channels(session_chList);
%         
%             % exclude EMG, reference channels
%             cp = initChanParams();
%             cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
%             sessionChannels = excludeChannels(cp, sessionChannels);

%             session_regionList = getRegionsfromChannelDB(sessionChannels);
%             num_sessionRegions = length(session_regionList);
%         
%             totalSessionChannels = length(sessionChannels);

            sessionRegionList = region_phaseAmp_metadata.regionList;
            num_sessionRegions = length(sessionRegionList);

            for iRegion = 1 : num_sessionRegions
            
%                 cp = initChanParams();
%                 cp.locationName = regionList{iRegion};            
%                 chList = extractChannels(cp, sessionChannels);
%                 if isempty(chList); continue; end
%                 regionChannels = sessionChannels(chList);
%                 numCh = length(regionChannels);

                % find index of current region in list of all regions across
                % sessions
                currentRegion = sessionRegionList{iRegion};
                all_regionIdx = find(strcmpi(currentRegion, allRegionList));
           
                mean_regionMI(all_regionIdx, :, :, :, :) = mean_regionMI(all_regionIdx, :, :, :, :) + mean_MI(iRegion, :, :, :, :);
                num_regionSessions(all_regionIdx) = num_regionSessions(all_regionIdx) + 1;
                num_channels_per_regionSession(iSession, all_regionIdx) = region_phaseAmp_metadata.num_ch_per_region(iRegion);
                
            end    % end for iRegion...
            
        end    % end for iSession...
        
        eventtWin = region_phaseAmp_metadata.eventtWin;
        tmin = eventtWin(1);
        analysisWin = region_phaseAmp_metadata.analysisWin;
        stepSize = region_phaseAmp_metadata.stepSize;

        figProps.n = numEventTypes;
        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
        numPages = 0;
        numRegionPlots = 0;
        for i_allRegions = 1 : num_allRegions
            mean_regionMI(i_allRegions, :, :, :, :) = mean_regionMI(i_allRegions, :, :, :, :) / num_regionSessions(i_allRegions);
        
            if makePlots
                numRegionPlots = numRegionPlots + 1;
                rowNum = rem(numRegionPlots, regions_per_page);
                
                if rowNum == 1
                    [h_fig, h_axes] = createFigPanels5(figProps);
                    page_regionList = allRegionList{i_allRegions};
                    page_numSessionList = num2str(num_regionSession(i_allRegions));
                    numPages = numPages + 1;
                else
                    page_regionList = [page_regionList ', ' allRegionList{i_allRegions}];
                    page_numSessionList = [page_numSessionList ', ' num2str(num_regionSession(i_allRegions))];
                end
                if rowNum == 0; rowNum = regions_per_page; end
                
                for iEventType = 1 : numEventTypes
                    
                    toPlot = log10(squeeze(mean_regionMI(i_allRegions, iEventType, :, phase_freq_idx, :))');
                    axes(h_axes(rowNum, iEventType));
                    imagesc(low_freqs, high_freqs, toPlot);
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
                    h_figAxes = createFigAxes(h_fig);
                    axes(h_figAxes);
                    
                    textStr{1} = ['phase-amplitude coupling, phase frequency = ' num2str(phase_freq)];
                    textStr{2} = ['Trial type: ' phaseAmp_metadata.trialType];
                    textStr{3} = page_regionList;
                    textStr{4} = ['Number of sessions for each region: ' page_numSessionList];
                    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                    
                    if numPages == 1
                        export_fig(fig_saveName,'-pdf','-q101','-painters');
                    else
                        export_fig(fig_saveName,'-pdf','-q101','-painters','-append');
                    end
                    close(h_fig);

                end    % if (rowNum == regions_per_page || numRegionPlots == num_allRegions)

            end    % if makePlots...
        
        end    % for i_allRegions = 1 : num_allRegions
        
        save(region_mean_saveName, 'mean_regionMI', 'num_regionSessions', 'all_region_phaseAmp_metadata');

    end

end   % for i_chDB = ...

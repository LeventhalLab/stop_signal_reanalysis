% script_plotFries_sfc

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
sfcRoot = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Fries_sfc';

[chDB_list, chDB_fnames, spikeDB_list, spikeDB_fnames] = get_chStructs_for_analysis();
numEventTypes = 5;

colorLim = [0 0.4];
                                          
for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant spike and channel DBs, if necessary
%     if ~exist(chDB_list{i_chDB}, 'var')
%         chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
%         disp(['loading ' chDB_file]);
%         load( chDB_file );
%     end
    if ~exist(spikeDB_list{i_chDB}, 'var')
        spikeDB_file = fullfile(chDB_directory, spikeDB_fnames{i_chDB});
        disp(['loading ' spikeDB_file]);
        load( spikeDB_file );
    end
% NOT SURE THE FULL CHANNEL/SPIKE DB'S ARE NEEDED

    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
        spikeDB_info = whos( [spikeDB_list{i_chDB}(1:3) 'spike'] );
    else
        implantID = chDB_list{i_chDB}(1:5);
        chDB_info = whos( [implantID 'Ch*'] );
        spikeDB_info = whos( implantID );
    end
%     channels = eval( chDB_info.name );
    spike    = eval( spikeDB_info.name );
    
    subject_sfcDir = fullfile(sfcRoot, [implantID, '_Fries_sfc']);
    cd(subject_sfcDir);
    
    sfcFiles = dir([implantID '*.mat']);
    numFiles = length(sfcFiles);
    
    for iFile = 1 : numFiles
        
        load(sfcFiles(iFile).name);
        numEventTypes = length(Fsfc_metadata.eventList);
        numUnits = length(spName);
        numSteps = floor(range(Fsfc_metadata.eventWin) / Fsfc_metadata.stepSize);
        
        cp = initChanParams();
        cp.channelName = spName{1};
        spList = extractChannels(cp, spike);
        sp = spike{spList};
        
        t1 = Fsfc_metadata.eventWin(1) + Fsfc_metadata.stepSize*0.5;
        t2 = t1 + (numSteps-1) * Fsfc_metadata.stepSize;
        t = linspace(t1,t2,numSteps);
        f = Fsfc_metadata.f;

        figProps.width  = 11 * 2.54;
        figProps.height = 8.5 * 2.54;

        figProps.m = 2; figProps.n = numEventTypes;   % number of rows and columns, respectively

        sideMargins = 2; botMargin = 2.54;

        figProps.colSpacing = ones(1, figProps.n) * 0.5;
        figProps.rowSpacing = ones(1, figProps.m) * 0.5;
        figProps.topMargin  = 5;

        figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                                     2 * sideMargins - ...
                                                     sum(figProps.colSpacing)) / figProps.n;
        figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                      figProps.topMargin - botMargin - ...
                                                      sum(figProps.rowSpacing)) / figProps.m;
                                                  
        saveName = [implantID, '_', sp.session, '_Fsfc.pdf'];
        saveName = fullfile(subject_sfcDir, saveName);
                                                  
        for iUnit = 1 : numUnits
            
            cp = initChanParams();
            cp.channelName = spName{iUnit};
            spList = extractChannels(cp, spike);
            sp = spike{spList};
                                          
            [h_fig, h_axes] = createFigPanels5(figProps);

            maxFreq = 0;
            for iEventType = 1 : numEventTypes

                im_to_plot = squeeze(sfc(iUnit, iEventType, :, :))';
                axes(h_axes(1, iEventType));
                imagesc(t,f,im_to_plot);
                set(gca,'ydir','normal','clim',colorLim);
                title(Fsfc_metadata.eventList{iEventType});
                if iEventType == 1
                    ylabel('frequency (Hz)');
                    colorbar
                else
                    set(gca,'yticklabel',[]);
                end
                
                axes(h_axes(2, iEventType));
                spikeFreq = (squeeze(numSpikes(iUnit, iEventType, :)) / Fsfc_metadata.numTrials) / Fsfc_metadata.spikeWin;
                plot(t, spikeFreq);
                maxFreq = max(maxFreq, max(spikeFreq));
                if iEventType == 1
                    ylabel('spike frequency');
                else
                    set(gca,'yticklabel',[]);
                end
            end
            for iEventType = 1 : numEventTypes
                axes(h_axes(2, iEventType));
                set(gca,'ylim',[0 maxFreq]);
            end
            
            h_figAxes = createFigAxes(h_fig);
            axes(h_figAxes);
            totSpikes = num2str(length(sp.timestamps.spike));
            textStr{1} = [spName{iUnit} ', ' sp.location.name '; ' totSpikes ' total spikes; ' num2str(Fsfc_metadata.numTrials) ' trials'];
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
            
            if iUnit == 1
                export_fig(saveName,'-pdf','-q101','-painters')
            else
                export_fig(saveName,'-pdf','-q101','-painters', '-append')
            end
            close(h_fig);
        
        end    % for iUnit...
        

    end    % end for iFile...
    
end
                
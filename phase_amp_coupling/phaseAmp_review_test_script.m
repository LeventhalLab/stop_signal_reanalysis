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
            

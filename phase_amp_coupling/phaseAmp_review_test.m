    %% 
%         mean_mrl_byRegio÷n = squeeze(nanmean(mean_mrl_acrossSessions, 1));
        mean_mrl_z_byRegion = squeeze(nanmean(mean_mrl_z_acrossSessions, 1));

        numPages = 0;
        numRegionPlots = 0;
        h_fig = cell(length(var_to_plot), length(plotTypes)); h_axes = cell(length(var_to_plot), length(plotTypes));
%%
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
                page_regionList = [page_regionList ', ' ch.name];
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

                            if (rowNum == regions_per_page || numChPlots == totalSessionChannels)
                                h_figAxes = createFigAxes(h_fig{iVar, iPlotType}(iFreq));
                                axes(h_figAxes);

                                textStr{2} = ['Trial type: ' region_phaseAmp_metadata.trialType];
                                textStr{3} = page_regionList;
                                textStr{4} = ['Number of sessions in average: ' numSessions_perRegionList];
                                text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
                                if numPages == 1
                                    export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters');
                                else
                                    export_fig(fig_saveName{iVar, iPlotType}{iFreq},'-pdf','-q101','-painters','-append');
                                end
                                close(h_fig{iVar, iPlotType}(iFreq));
                            end

                        end    % for iEventType...
                    end    % for iFreq...
                end    % for iPlotType...

            end    % for iVar...

        end    % for iRegion...
            

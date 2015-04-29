% script to plot power spectra for each region

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
powerSpectraDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectra';
lfp_root          = ['/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs'];
sample_lfp_name = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs/IM164_LFPs/IM164_20091112_13-56-32.hsdf';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);
flength = 1025;    % just to preallocate matrices

numRows = 4;
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;
sideMargins = 2; botMargin = 2.54;
figProps.topMargin  = 4;
figProps.n = numRegions;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.colSpacing = ones(1, figProps.n) * 0.5;
figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                             2 * sideMargins - ...
                                             sum(figProps.colSpacing)) / figProps.n;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

ylim = [-15 0];
max_f_to_plot = [10,100];
for i_chDB = 10:11%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
%     if ~exist(chDB_list{i_chDB}, 'var')
%         chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
%         disp(['loading ' chDB_file]);
%         load( chDB_file );
%     end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    
    subject_powerSpectraDir = fullfile(powerSpectraDir, [implantID '_ps']);
    if ~exist(subject_powerSpectraDir,'dir')
        continue
    end
    PDFname = fullfile(subject_powerSpectraDir, [implantID '_ps_plots.pdf']);
    
    cd(subject_powerSpectraDir);
    fileList = dir('ps_*');
    numSessions = length(fileList);
    
%     if i_chDB < 5
%         chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
%     else
%         chDB_info = whos( [implantID 'Ch*'] );
%     end
%     channels = eval( chDB_info.name );
%     
%     sessionList = getSessionsfromChannelDB( channels );
%     numSessions = length( sessionList );
    
    mean_session_psd = NaN(numSessions, numRegions, flength);
    numPages = 0;
    for iSession = 1 : numSessions
        figProps.m = numRows; figProps.n = numRegions;
        figProps.rowSpacing = ones(1, figProps.m) * 0.5;
        figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                                  figProps.topMargin - botMargin - ...
                                                  sum(figProps.rowSpacing)) / figProps.m;
        fname = fileList(iSession).name;
        load(fname);
        
        plotRow = rem(iSession,figProps.m);
        if plotRow == 1
            [h_fig, h_axes] = createFigPanels5(figProps);
            page_sessionList = ps_metadata.session;
            numPages = numPages + 1;
        else
            page_sessionList = [page_sessionList ', ' ps_metadata.session];
        end
        if plotRow == 0
            plotRow = figProps.m;
        end
        
        numCh = length(ps_metadata.region);
        for iRegion = 1 : numRegions
            
            axes(h_axes(plotRow,iRegion))
            
            numRegionChannels = 0;regionChList = [];
            for iCh = 1 : numCh
                if strcmpi(ps_metadata.region{iCh},ROI_list{iRegion})
                    numRegionChannels = numRegionChannels + 1;
                    regionChList = [regionChList,iCh];
                end
            end
            if numRegionChannels > 0
                session_ch_means = NaN(numRegionChannels, flength);

                for iRegionCh = 1 : numRegionChannels
                    session_ch_means(iRegionCh,:) = mean(squeeze(pxx(regionChList(iRegionCh),:,:)),2)';
                    plot(f(f<max_f_to_plot(2)), log(session_ch_means(iRegionCh,f<max_f_to_plot(2))),'k')
                    hold on
                    
                end
                mean_session_psd(iSession, iRegion, :) = nanmean(session_ch_means,1);
                plot(f(f<max_f_to_plot(2)), log(squeeze(mean_session_psd(iSession,iRegion,f<max_f_to_plot(2)))),'r')
            end
            set(gca,'ylim',ylim);
            if plotRow == 1
                title(ROI_list{iRegion});
            end
            if plotRow < figProps.m
                set(gca,'xtick',[],'xticklabel',[]);
            else
                set(gca,'xtick',[0:20:max_f_to_plot(2)]);
            end
            
            if iRegion == 1
                ylabel(ps_metadata.session);
            else
            end
            
        end   % for iRegion...
        
        if plotRow == figProps.m || iSession == numSessions
            h_figAxes = createFigAxes(h_fig);
            axes(h_figAxes);
            
            textStr{1} = [implantID ', power spectra, averages within sessions'];
            text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

            if numPages == 1
                export_fig(PDFname, '-pdf', '-q101', '-painters','-nocrop');
            else
                export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
            end
            close(h_fig);
        end
            
    end    % for iSession...
    
    figProps.m = 2; figProps.n = numRegions;
    figProps.rowSpacing = ones(1, figProps.m) * 0.5;
    figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;
    [h_fig, h_axes] = createFigPanels5(figProps);
    h_figAxes = createFigAxes(h_fig);
    for ii = 1 : 2
        for iRegion = 1 : numRegions
            axes(h_axes(ii,iRegion));
            for iSession = 1 : numSessions
                plot(f(f<max_f_to_plot(ii)), log(squeeze(mean_session_psd(iSession,iRegion,f<max_f_to_plot(ii)))),'k')
                hold on
            end
            plot(f(f<max_f_to_plot(ii)), log(squeeze(nanmean(mean_session_psd(:,iRegion,f<max_f_to_plot(ii)), 1))),'r')
            set(gca,'ylim',ylim);
        end
        if ii == 1
            title(ROI_list{iRegion});
        end
    end
    
    axes(h_figAxes);
    textStr{1} = [implantID ', power spectra, averages across sessions'];
    export_fig(PDFname, '-pdf', '-q101', '-painters', '-append','-nocrop');
    text('units','centimeters','position',[3, 8*2.54], 'string',textStr);
    close(h_fig); 
    
end
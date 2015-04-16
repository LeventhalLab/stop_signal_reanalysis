% script_plot_STOPsuccvfail_z_across_subjects

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

regions_per_page = 5;
colorLim = [0 1];
vecDiff_colorLim = [0 1];
z_yLim = [-1 2];


sideMargins = 2; botMargin = 2.54;


figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

numSubs = 4;
ROI_list = {'EEG','Str','GP','STN','SNr'};
ROI_idx = NaN(numSubs, length(ROI_list));
numROI = zeros(1, length(ROI_list));

testFile = fullfile(phase_stop_succvfail_directory, 'IM166_stopPhases', 'IM166_vecDiff_z_stopSuccvFail.mat');
load(testFile);
numEventTypes = length(STOPmrv_z_acrossSessions_metadata(1).eventList);
numFreq = size(STOPmrv_z_acrossSessions_metadata(1).freqBands,1);
z = zeros(numSubs,length(ROI_list), numEventTypes, size(region_z,3), size(region_z, 4));

regions_per_page = 5;

figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = regions_per_page;    % number of rows
figProps.n = numEventTypes;

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

t = linspace(-1,1,size(region_z,4));
f = mean(STOPmrv_z_acrossSessions_metadata(1).freqBands,2);
phase_targetFreq = 2.5;
for i_chDB = 1:numSubs%length(chDB_list)
    i_chDB
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        disp([subject_stopPhaseDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    
    vecDiff_z_acrossSessions_mat_saveName = [implantID '_vecDiff_z_stopSuccvFail.mat'];
    vecDiff_z_acrossSessions_mat_saveName = fullfile(subject_stopPhaseDir, vecDiff_z_acrossSessions_mat_saveName);
    
    load(vecDiff_z_acrossSessions_mat_saveName);
    
    for iROI = 1 : length(ROI_list)
        temp = find(strcmpi(ROI_list{iROI}, allRegionList));
        if ~isempty(temp)
            ROI_idx(i_chDB, iROI) = temp;
            
            z(i_chDB,iROI,:,:,:) = squeeze(z(i_chDB,iROI,:,:,:)) + ...
                                   squeeze(region_z(ROI_idx(i_chDB, iROI),:,:,:));
            numROI(iROI) = numROI(iROI) + 1;
        end
    end
    
end

zphase_acrossSubjects_pdfname = fullfile(phase_stop_succvfail_directory, 'z_stop_succvfail_across_subjects_singlefreq.pdf');

[h_fig, h_axes] = createFigPanels5(figProps);

phase_freqIdx = find(abs(f - phase_targetFreq) == min(abs(f - phase_targetFreq)));
mean_z = squeeze(nanmean(z,1));
for iROI = 1 : length(ROI_list)

%     z(iROI,:,:,:) = z(iROI,:,:,:) / numROI(iROI);

    for iEventType = 1 : numEventTypes
        axes(h_axes(iROI, iEventType));

        z_toPlot = squeeze(mean_z(iROI,iEventType, phase_freqIdx, :));
        plot(t,z_toPlot,'color','r','linewidth',2);
        hold on
        for iSub = 1 : numSubs
            if ~isnan(ROI_idx(iSub, iROI))
                toPlot = squeeze(z(iSub, iROI, iEventType, phase_freqIdx, :));
                plot(t, toPlot, 'color','r','linestyle','--');
            end
        end
        line([-1,1],[0,0],'color','k','linestyle','--')
        line([0,0i],z_yLim,'color','k','linestyle','--')
        set(gca,'ydir','normal','ylim', z_yLim);
        if iROI == 1
            title(STOPmrv_z_acrossSessions_metadata(1).eventList{iEventType});
        end
        if iROI < regions_per_page
            set(gca,'xticklabel',[]);
        else
            xlabel('time (s)')
        end
        if iEventType > 1
            set(gca,'yticklabel',[]);
        else
            ylabel('z-score')
        end
    end
end

export_fig(zphase_acrossSubjects_pdfname,'-pdf','-q101','-painters');
    
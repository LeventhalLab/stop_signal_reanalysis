% script_plotRT_powerCorrelations

bitOrder = 'b';
colorLim = [-0.4, 0.4];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT correlation plots';

ROI_list = {'eegorb','cpu','gp','stn','snr'};                               % regions of interest
EOI_list = {'cueon','nosecenterin','tone','nosecenterout','nosesidein'};    % events of interest

%     
[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_powerRT_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations/IM164_powerRTcorr/D2220091005/power_RTcorr_D2220091005T02.mat';

load(sample_powerRT_file);
eventList = powerRTcorr_metadata.eventList;
numEvents = length(eventList);
if numEvents == 6    % working around having included the same event twice in the power-RT calculations
    if strcmpi(eventList{5},eventList{6})
        eventList = eventList(1:5);
        numEvents = 5;
    end
end
numSamps = size(powerRTcorr, 3);
numFreqs = size(powerRTcorr, 2);
twin = powerRTcorr_metadata.twin;
t = linspace(twin(1), twin(2), numSamps);
f = powerRTcorr_metadata.freqs;

xticks = [twin(1),0,twin(2)];
yticks = 20:20:100;%min(f):20:max(f);
    
% try making a separate pdf for each session
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page; figProps.n = length(EOI_list);   % number of rows and columns, respectively

sideMargins = 1.5; botMargin = 2.54;

figProps.colSpacing = ones(1, figProps.n) * 0.5;
figProps.rowSpacing = ones(1, figProps.m) * 0.5;
figProps.topMargin  = 5;

figProps.panelWidth = ones(1, figProps.n) * (figProps.width - ...
                                             2 * sideMargins - ...
                                             sum(figProps.colSpacing)) / figProps.n;
figProps.panelHeight = ones(1, figProps.m) * (figProps.height - ...
                                              figProps.topMargin - botMargin - ...
                                              sum(figProps.rowSpacing)) / figProps.m;

numSubs = 3;
RTcorr_subjects = zeros(numSubs, length(ROI_list), length(eventList), numFreqs, numSamps);
numROI = zeros(1, length(ROI_list));
ROI_idx = zeros(1, length(ROI_list));
for i_chDB = 1 : numSubs%length(chDB_list)
    
%     % first, load the relevant channel DBs, if necessary
%     if ~exist(chDB_list{i_chDB}, 'var')
%         chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
%         disp(['loading ' chDB_file]);
%         load( chDB_file );
%     end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    
    subject_powerRTcorr_directory = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr']);
    if ~exist(subject_powerRTcorr_directory, 'dir')
        disp([subject_powerRTcorr_directory ' not found.']);
        continue;
    end

    subject_RTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_plots']);
    if ~exist(subject_RTcorr_plots_directory, 'dir')
        mkdir(subject_RTcorr_plots_directory);
    end
    
    regionSummaryMatName = [implantID '_powerRTcorr_across_sessions.mat'];
    regionSummaryMatName = fullfile(subject_powerRTcorr_directory, regionSummaryMatName);
    if ~exist(regionSummaryMatName, 'file');continue;end
    
    load(regionSummaryMatName);
    
    % find indices of regions we care about
    for iROI = 1 : length(ROI_list)
        temp = find(strcmpi(ROI_list{iROI}, region_power_RTcorr_metadata.regionList));
        if isempty(temp)
            ROI_idx(iROI) = NaN;
        else
            ROI_idx(iROI) = temp;
            numROI(iROI)  = numROI(iROI) + 1;
        end
    end
    
    % find indices of events we care about
    EOI_idx = zeros(1, length(EOI_list));
    for iEOI = 1 : length(EOI_list)
        temp = find(strcmpi(EOI_list{iEOI}, region_power_RTcorr_metadata.eventList));
        if isempty(temp)
            EOI_idx(iEOI) = NaN;
        else
            EOI_idx(iEOI) = temp;
        end
    end
    
    for iROI = 1 : length(ROI_list)
        for iEOI = 1 : length(EOI_list)
            if ~isnan(ROI_idx(iROI))
                RTcorr_subjects(i_chDB, iROI, iEOI, :, :) = squeeze(RTcorr_subjects(i_chDB, iROI, iEOI, :, :)) + ...
                                                                   squeeze(mean_RTcorr(ROI_idx(iROI), EOI_idx(iEOI), :, :));
            end
        end
    end
    
end

PDFname = fullfile(RTcorr_plots_directory, 'powerRTcorr_across_subjects.pdf');
[h_fig, h_axes] = createFigPanels5(figProps);
RTcorr_across_subjects = zeros(length(ROI_list), length(eventList), numFreqs, numSamps);
regionString = '';
numSubjectsString = '';
for iROI = 1 : length(ROI_list)
    
    regionString = [regionString ', ' ROI_list{iROI}];
    numSubjectsString = [numSubjectsString ', ' num2str(numROI(iROI))];
    for iEOI = 1 : length(EOI_list)
        temp = sum(squeeze(RTcorr_subjects(:, iROI, iEOI, :, :)), 1) / numROI(iROI);
        
        RTcorr_across_subjects(iROI, iEOI, : ,:) = squeeze(temp);

        axes(h_axes(iROI, iEOI));
        
        toPlot = squeeze(RTcorr_across_subjects(iROI, iEOI, :, :));
        imagesc(t, f, toPlot);
        set(gca,'ydir','normal');
        set(gca,'clim',colorLim);
        if iROI == 1
            title(eventList{EOI_idx(iEOI)});
        end
        
        if iROI < figProps.m
            set(gca,'xticklabel',{});
        else
            xlabel('time (s)')
            set(gca,'xtick',xticks);
        end
        if iEOI > 1
            set(gca,'yticklabel',{});
        else
            ylabel('frequency (Hz)')
            set(gca,'ytick',yticks);
        end

    end    % for iEOI
    
end    % for iROI

colorbar

h_figAxes = createFigAxes(h_fig);
axes(h_figAxes);

textStr{1} = ['Regions: ' regionString];
textStr{2} = ['number of subjects per region: ' numSubjectsString];
textStr{5} = ['color limits: ' num2str(colorLim(1)) ' to ' num2str(colorLim(2))];
text('units','centimeters','position',[3, 8*2.54], 'string',textStr);

export_fig(PDFname, '-pdf', '-q101', '-painters');

% close(h_fig);

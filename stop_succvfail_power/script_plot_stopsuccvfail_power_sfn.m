% script_plotRT_powerCorrelations

bitOrder = 'b';
colorLim = [-3, 3];
xLim = [-1, 1];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT correlation plots';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

ROI_list = {'eeg','str','gp','stn','snr'};                               % regions of interest
EOI_list = {'cueon','nosecenterin','tone','whiteNoise','noseCenterOut'};    % events of interest

%     
[chDB_list, chDB_fnames] = get_chStructs_for_analysis;
channels_per_page = 5;

% load a sample file
sample_succvfail_file = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power/IM166_stopPower/IM166_power_stopSuccvFail_z.mat';


load(sample_succvfail_file);
eventList = STOPpower_z_acrossSessions_metadata.eventList;
numEvents = length(eventList);
if numEvents == 6    % working around having included the same event twice in the power-RT calculations
    if strcmpi(eventList{5},eventList{6})
        eventList = eventList(1:5);
        numEvents = 5;
    end
end
numSamps = size(region_z, 4);
numFreqs = size(region_z, 3);
twin = STOPpower_z_acrossSessions_metadata.twin;
t = linspace(twin(1), twin(2), numSamps);
f = STOPpower_z_acrossSessions_metadata.f;

xticks = [xLim(1),0,xLim(2)];
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

numSubs = 4;
powerSuccFail_subjects = zeros(numSubs, length(ROI_list), length(eventList), numFreqs, numSamps);
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
    
    subject_stopsuccvfail_directory = fullfile(power_stop_succvfail_directory, [implantID '_stopPower']);
    if ~exist(subject_stopsuccvfail_directory, 'dir')
        disp([subject_stopsuccvfail_directory ' not found.']);
        continue;
    end

%     subject_RTcorr_plots_directory = fullfile(power_stop_succvfail_directory, [implantID '_RTcorr_plots']);
%     if ~exist(subject_RTcorr_plots_directory, 'dir')
%         mkdir(subject_RTcorr_plots_directory);
%     end
%     
    regionSummaryMatName = [implantID '_power_stopSuccvFail_z.mat'];
    regionSummaryMatName = fullfile(subject_stopsuccvfail_directory, regionSummaryMatName);
    if ~exist(regionSummaryMatName, 'file');continue;end
    
    load(regionSummaryMatName);
    
    % find indices of regions we care about
    for iROI = 1 : length(ROI_list)
        temp = find(strcmpi(ROI_list{iROI}, STOPpower_z_acrossSessions_metadata.regions));
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
        temp = find(strcmpi(EOI_list{iEOI}, STOPpower_z_acrossSessions_metadata.eventList));
        if isempty(temp)
            EOI_idx(iEOI) = NaN;
        else
            EOI_idx(iEOI) = temp;
        end
    end
    
    for iROI = 1 : length(ROI_list)
        for iEOI = 1 : length(EOI_list)
            if ~isnan(ROI_idx(iROI))
                powerSuccFail_subjects(i_chDB, iROI, iEOI, :, :) = squeeze(powerSuccFail_subjects(i_chDB, iROI, iEOI, :, :)) + ...
                                                                   squeeze(region_z(ROI_idx(iROI), EOI_idx(iEOI), :, :));
            end
        end
    end
    
end

PDFname = fullfile(power_stop_succvfail_directory, 'power_succvfail_across_subjects_4.pdf');
[h_fig, h_axes] = createFigPanels5(figProps);
z_across_subjects = zeros(length(ROI_list), length(eventList), numFreqs, numSamps);
regionString = '';
numSubjectsString = '';
for iROI = 1 : length(ROI_list)
    
    regionString = [regionString ', ' ROI_list{iROI}];
    numSubjectsString = [numSubjectsString ', ' num2str(numROI(iROI))];
    for iEOI = 1 : length(EOI_list)
        temp = sum(squeeze(powerSuccFail_subjects(:, iROI, iEOI, :, :)), 1) / numROI(iROI);
        
        z_across_subjects(iROI, iEOI, : ,:) = squeeze(temp);

        axes(h_axes(iROI, iEOI));
        
        toPlot = squeeze(z_across_subjects(iROI, iEOI, :, :));
        imagesc(t, f, toPlot);
        set(gca,'ydir','normal','xlim',xLim);
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

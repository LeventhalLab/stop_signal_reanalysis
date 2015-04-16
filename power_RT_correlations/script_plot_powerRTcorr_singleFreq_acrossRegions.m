bitOrder = 'b';
yLim = [-1 1];

chDB_directory         = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
powerRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations';
RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT correlation plots';

regions_of_interest = {'eegorb','Cpu','GP','STN','SNr'};

regions_per_page = 5;

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
yticks = [yLim(1), 0, yLim(2)];

% try plotting all animal averages on single graphs, including the averages
figProps.width  = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.m = channels_per_page; figProps.n = numEvents;   % number of rows and columns, respectively

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
                                          
[h_fig, h_axes] = createFigPanels5(figProps);
                                          
targetFreq = 19.5;
freqStr = [num2str(round(targetFreq)) 'Hz'];
                                          
for i_chDB = 1 : length(chDB_list)
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_powerRTcorr_directory = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr']);
    subject_powerRTcorr_plots_directory = fullfile(RTcorr_plots_directory, [implantID '_RTcorr_plots']);
    
    regionSummaryMatName = [implantID '_powerRTcorr_across_sessions.mat'];
    regionSummaryMatName = fullfile(subject_powerRTcorr_directory, regionSummaryMatName);
    if ~exist(regionSummaryMatName,'file');continue;end
    
    load(regionSummaryMatName);
    
    PDFname = fullfile(subject_RTcorr_plots_directory, [implantID '_powerRTcorr_' freqStr '.pdf']);
    
    f = region_power_RTcorr_metadata.f;
    
    freqIdx = find(abs(f - targetFreq) == min(abs(f - targetFreq)));
    
    
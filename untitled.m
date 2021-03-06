% script_plotMI

% script to plot the modulation index heat maps to look for phase-amplitude
% coupling

bitOrder = 'b';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
comod_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase-amplitude comodugrams';
comod_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase-amplitude comod plots';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_comoddir = fullfile(comod_directory, [implantID '_phase-amp_comods']);
    if ~exist(subject_comoddir, 'dir')
        disp([subject_comoddir ' not found.']);
        continue;
    end
    
    subject_comodPlotsDir = fullfile(comod_plots_directory, [implantID '_comodPlots']);
    if ~exist(subject_comodPlotsDir, 'dir')
        mkdir(subject_comodPlotsDir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        numCh = length(sessionChannels);
        
        comod_sessionDir = fullfile(subject_comoddir, sessionList{iSession});
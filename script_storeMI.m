% script_storeMI

bitOrder = 'b';
low_freq_range  = [0 25];
high_freq_range = [13 101];

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
comod_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase-amplitude comodugrams';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

for i_chDB = 2 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir, 'dir')
        disp([subject_hilbertDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    subject_comoddir = fullfile(comod_directory, [implantID '_phase-amp_comods']);
    if ~exist(subject_comoddir, 'dir')
        mkdir(subject_comoddir);
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
        
        hilbert_sessionDir = fullfile(subject_hilbertDir, sessionList{iSession});
        comod_sessionDir   = fullfile(subject_comoddir, sessionList{iSession});
        if ~exist(comod_sessionDir, 'dir')
            mkdir(comod_sessionDir);
        end
        
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
        
        if ~exist(metadata_filename, 'file')
            continue
%             error([metadata_filename ' could not be found.']);
        end
        load(metadata_filename);
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            comod_name = [ch.name '_phase_amp_comodugram.mat'];
            comod_name = fullfile(comod_sessionDir, comod_name);
            
            if exist(comod_name, 'file')
                continue;
            end
        
            [MI, low_freqs, high_freqs] = phase_amp_comodugram(sessionChannels{iCh}, ...
                                                               'lowfreqrange', low_freq_range, ...
                                                               'highfreqrange', high_freq_range);
                                                           
            save(comod_name, 'MI', 'low_freqs', 'high_freqs', 'metadata');
            
        end
        
    end
    
end
            
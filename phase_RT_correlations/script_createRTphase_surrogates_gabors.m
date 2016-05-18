% script_createRTphase_surrogates_gabors
%
% script to create surrogate mrl distributions for analyzing relationship
% between phase and RT
%

bitOrder = 'b';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
gabor_directory =  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations_gabors';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

for i_chDB = 1 : length(chDB_list)
    
    % first,  the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['ing ' chDB_file]);
        ( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_gaborDir = fullfile(gabor_directory, [implantID '_gabor']);
    if ~exist(subject_gaborDir, 'dir')
        disp([subject_gaborDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    % exclude EMG, reference channels
    cp = initChanParams();
    cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
    channels = excludeChannels(cp, channels);
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iSession = 1 : numSessions
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        numCh = length(sessionChannels);
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            createSurrogate_phase_RT_dist(ch, ...
                'gabordirectory', gabor_directory, ...
                'phasertcorrdir', phaseRTcorr_directory);
            
        end
        
    end
    
end
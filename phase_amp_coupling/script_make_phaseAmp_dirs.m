trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;



for i_chDB = 3 : 3%length(chDB_list)
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
    for iTrialType = 1 : length(trialTypeList)
        trialType = trialTypeList{iTrialType};
        
        
        for iSession = 1 : numSessions
        
            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                mkdir(phaseAmp_sessionDir);
            end
            
        end
        
    end
    
end
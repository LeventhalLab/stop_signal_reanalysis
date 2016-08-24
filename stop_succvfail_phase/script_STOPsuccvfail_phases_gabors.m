% script_STOPsuccvfail_phases_gabors

% script to go through all STOP sessions and save continuous phase values
% for stop-success and stop-fail trials

% chDB_directory    = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% hilbert_1Hz_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
% phase_stop_succvfail_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase';
% power_stop_succvfail_directory = '\\141.214.45.212/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
gabor_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
% hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
% hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phase_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_phase_gabors';
power_stop_succvfail_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop_succvfail_power_gabors';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
numEventTypes = length(eventList);
twin = [-1 1];

trialTypeList = {'correctstop','failedstop'};

for i_chDB = 2 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
        
    subject_gaborDir = fullfile(gabor_directory, [implantID '_ps']);
%     subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
%     subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    if ~exist(subject_gaborDir, 'dir')
        disp([subject_gaborDir ' not found. Skipping ' implantID '...'])
        continue
    end
    
    subject_stopPhaseDir = fullfile(phase_stop_succvfail_directory, [implantID '_stopPhases']);
    if ~exist(subject_stopPhaseDir, 'dir')
        mkdir(subject_stopPhaseDir);
    end
    subject_stopPowerDir = fullfile(power_stop_succvfail_directory, [implantID '_stopPower']);
    if ~exist(subject_stopPowerDir, 'dir')
        mkdir(subject_stopPowerDir);
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    cp = initChanParams();
    cp.task = 3;
    chList = extractChannels(cp, channels);
    channels = channels(chList);
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
    for iSession = 1 : numSessions

        disp(sessionList{iSession});
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        cp.task = 3;    % only stop-signal sessions
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        cp = initChanParams();
        cp.tetrode = {'e2', 'e3','e02','e03'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numCh = length(sessionChannels);
        
        session_stopPhaseDir = fullfile(subject_stopPhaseDir, sessionList{iSession});
        if ~exist(session_stopPhaseDir, 'dir')
            mkdir(session_stopPhaseDir);
        end
        session_stopPowerDir = fullfile(subject_stopPowerDir, sessionList{iSession});
        if ~exist(session_stopPowerDir, 'dir')
            mkdir(session_stopPowerDir);
        end
        
        gabor_sessionDir = fullfile(subject_gaborDir, [sessionList{iSession} '_scalograms']);
%         hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
%         hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});
        if ~exist(gabor_sessionDir, 'dir')
            disp([gabor_sessionDir ' could not be found.']);
            continue;
        end

        for iCh = 1 : numCh
%             iCh
            
            ch = sessionChannels{iCh};
            
            saveName_phase = ['stopPhases_' ch.name '_gabor.mat'];
            saveName_phase = fullfile(session_stopPhaseDir, saveName_phase);
%             if exist(saveName_phase, 'file'); continue; end
            
            saveName_power = ['stopPower_' ch.name '_gabor.mat'];
            saveName_power = fullfile(session_stopPowerDir, saveName_power);
            if exist(saveName_power, 'file') && exist(saveName_phase, 'file'); continue; end
            
            sclgrms = cell(1,2);
            validGaborsFound = true;
            for iTrialType = 1 : 2
                trialType = trialTypeList{iTrialType};
                gabor_name = [ch.name '_' trialType '_scalograms.mat'];
                gabor_name = fullfile(gabor_sessionDir, gabor_name);
                if ~exist(gabor_name,'file');
                    validGaborsFound = false;
                    break;
                end

                sclgrms{iTrialType} = load(gabor_name);
            end
            
            if ~validGaborsFound; continue; end

            STOPmetadata.Fs = sclgrms{1}.scalogram_metadata.Fs;
            STOPmetadata.f = sclgrms{1}.scalogram_metadata.f;
            STOPmetadata.t = sclgrms{1}.scalogram_metadata.t;
            
            STOPmetadata.eventList = eventList;
            STOPmetadata.twin      = sclgrms{1}.scalogram_metadata.twin;
            
            STOPmetadata.chName = sclgrms{1}.scalogram_metadata.channel;
            
            failedSTOP_power = abs(sclgrms{2}.W) .^2;
            correctSTOP_power = abs(sclgrms{1}.W) .^2;

            failedSTOP_phase = angle(sclgrms{2}.W);
            correctSTOP_phase = angle(sclgrms{1}.W);
            
            save(saveName_phase, 'correctSTOP_phase', 'failedSTOP_phase', 'STOPmetadata');
            save(saveName_power, 'correctSTOP_power', 'failedSTOP_power', 'STOPmetadata');
            
        end    % for iCh...
        
    end    % for iSession...
    
end    % for i_chDB
                    

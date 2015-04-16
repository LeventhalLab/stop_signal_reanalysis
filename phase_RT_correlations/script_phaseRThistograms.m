% script_phaseRThistograms


chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
numEvents = length(eventList);
twin = [-1 1];

trialType = 'correctgo';

for i_chDB = 1 : length(chDB_list)
    
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
    
    subject_phaseRThistdir = fullfile(phaseRThist_directory, [implantID '_phaseRTcorr']);
    if ~exist(subject_phaseRThistdir, 'dir')
        mkdir(subject_phaseRThistdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length(sessionList);
    
%     RT = cell(1, numSessions);
%     freqPhase = cell(1, numSessions);
    for iSession = 1 : numSessions

        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        numCh = length(sessionChannels);

        hilbert_sessionDir   = fullfile(subject_hilbertDir, sessionList{iSession});
        if ~exist(hilbert_sessionDir, 'dir'); continue; end
        
        phaseRThist_sessionDir = fullfile(subject_phaseRThistdir, sessionList{iSession});
        if ~exist(phaseRThist_sessionDir, 'dir')
            mkdir(phaseRThist_sessionDir);
        end
        
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
        
        if ~exist(metadata_filename, 'file')
            error([metadata_filename ' could not be found.']);
        end
        load(metadata_filename);
        
        centerFreqs = mean(metadata.freqBands, 2);
        numFreqs = length(centerFreqs);
        samps_per_window = round(range(twin) * metadata.Fs);
        
        phaseRThist_metadata.freqs = centerFreqs;
        phaseRThist_metadata.eventList = eventList;
        phaseRThist_metadata.twin = twin;
        phaseRThist_metadata.Fs = metadata.Fs;
        phaseRThist_metadata.trialType = trialType;
        
        numSamps = round(range(twin) * metadata.Fs);
        
        [tempRT, ~, ~] = collect_RT_MT_by_rat(sessionChannels, trialType);
        RT = tempRT{1};
        
%         tic
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            hilbert_name = ['analytic_' sessionChannels{iCh}.name '.bin'];
            hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
            
            phase_RThist_name = ['phase_RT_hist_analysis_' sessionChannels{iCh}.name '.mat'];
            phase_RThist_name = fullfile(phaseRThist_sessionDir, phase_RThist_name);
            if exist(phase_RThist_name, 'file'); continue; end
            
            phaseRThist_metadata.channel = ch.name;
        
%             tic
            freqPhase = zeros(numEvents, numFreqs, length(RT), samps_per_window);
            for iEvent = 1 : numEvents
%                 iEvent
%                 tic
                for iFreq = 1 : numFreqs
                    
                    ansig = getAnalyticAroundEvent( sessionChannels{iCh}, ...
                                                    iFreq, ...
                                                    eventList{iEvent}, ...
                                                    trialType, ...
                                                    twin );
                                                
                    freqPhase(iEvent, iFreq, :, :) = angle(ansig);
                    
                    
                end
                
            end
%             toc
            save(phase_RThist_name, 'RT', 'freqPhase', 'phaseRThist_metadata');
            
        end
%         toc
        
        
    end
    
end
                    
            
            
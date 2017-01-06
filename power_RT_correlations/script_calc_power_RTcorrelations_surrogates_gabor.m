% script_calc_power_RTcorrelations_surrogates_gabor

% UPDATED 10-27-2014

% ALSO NEED TO COMPARE PHASE OF ONGOING OSCILLATIONS IN STOP-SUCCESS VS
% STOP-FAILURE TRIALS AND NOGO-SUCCESS VS NOGO-FAIL TRIALS

% chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
% gabor_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
% powerRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_RT_correlations_gabors';
% lfp_root          = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/high_cutoff_stop-signal LFPs';

chDB_directory    = '/Volumes/Tbolt_02/stop-signal reanalysis/stop-signal data structures';
gabor_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/trial_scalograms';
powerRTcorr_directory = '/Volumes/Tbolt_02/stop-signal reanalysis/power_RT_correlations_gabors';
lfp_root          = '/Volumes/Tbolt_02/stop-signal reanalysis/high_cutoff_stop-signal LFPs';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

% eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
% numEvents = length(eventList);
twin = [-1 1];

trialType = 'correctgo';
numSurrogates = 200;

for i_chDB = 1 : 4%length(chDB_list)

    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_gabor_dir = fullfile(gabor_directory, [implantID '_ps']);
    subject_lfp_dir = fullfile(lfp_root, [implantID '_HF_LFPs']);
    if ~exist(subject_gabor_dir, 'dir')
        disp([subject_gabor_dir ' not found. Skipping ' implantID '...'])
        continue
    end
    % get names of all the lfp files
    cd(subject_lfp_dir);
    lfpDirStruct = dir([implantID '*.hsdf']);
    lfpList = cell(1,length(lfpDirStruct));
    for ii = 1 : length(lfpDirStruct)
        lfpList{ii} = lfpDirStruct(ii).name;
    end
    
    subject_powerRTcorrdir = fullfile(powerRTcorr_directory, [implantID '_powerRTcorr_gabor']);
    if ~exist(subject_powerRTcorrdir, 'dir')
        mkdir(subject_powerRTcorrdir);
    end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    [RT, ~, sessionList] = collect_RT_MT_by_rat(channels, trialType);
    numSessions = length( sessionList );
    
    % establish RT quantiles for phase analysis
%     allRT = RT{1};
%     if numSessions > 1
%         for iSession = 2 : numSessions
%             allRT = [allRT, RT{iSession}];
%         end
%     end
%     allRT = sort(allRT);
%     numRT = length(allRT);
    
    powerRTcorr = cell(1, numSessions);    % cell array to store correlations between continuous narrow-band power and RT

    for iSession = 1 : numSessions
        
        if iSession > 33 && iSession < 40;continue;end

        cp = initChanParams();
        cp.session = sessionList{iSession};
        
        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);
        
        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        
        cp = initChanParams();
        cp.tetrode = {'e2', 'e3'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numCh = length(sessionChannels);
        if numCh == 0; continue;end
        
        sessionDate = sessionChannels{1}.date;
        sessionDate = strrep(sessionDate,'-','');
        
        lfp_name_idx = strfind(lfpList,sessionDate);
        valid_lfp_name = false(1,length(lfpList));
        for ii = 1 : length(lfpList)
            valid_lfp_name(ii) = ~isempty(lfp_name_idx{ii});
        end
        lfpDuration = 0;
        if sum(valid_lfp_name) > 1
            disp(['more than one lfp file for ' sessionDate]);
        elseif sum(valid_lfp_name) == 0
            disp(['no lfp file for ' sessionDate]);
        else
            full_lfpName = fullfile(subject_lfp_dir, lfpList{valid_lfp_name});
            lfpDuration = getHSDlength( 'filename', full_lfpName );
        end

        gabor_sessionDir   = fullfile(subject_gabor_dir, [sessionList{iSession} '_scalograms']);
        if ~exist(gabor_sessionDir, 'dir'); continue; end
        
        powerRTcorr_sessionDir = fullfile(subject_powerRTcorrdir, sessionList{iSession});
        if ~exist(powerRTcorr_sessionDir, 'dir')
            mkdir(powerRTcorr_sessionDir);
        end
        
        powerRTcorr_surr_metadata.trialType = trialType;
        
        trialEventParams = getTrialEventParams('correctgo');
        correctGOtrials = extractTrials2(sessionChannels{1}.trials,trialEventParams);
        current_RT = RT{iSession};
        
        gabor_name = [sessionChannels{1}.name '_' trialType '_scalograms.mat'];
        gabor_name = fullfile(gabor_sessionDir, gabor_name);
        if ~exist(gabor_name,'file')
            disp([gabor_name ' not found.'])
            continue;
        end
        load(gabor_name);
        
        if correctGOtrials(1).timestamps.cueOn < (1-scalogram_metadata.twin(1))
            current_RT = current_RT(2:end);
        end
        if lfpDuration > 0
            if (lfpDuration - correctGOtrials(end).timestamps.noseSideOut) < (1+scalogram_metadata.twin(2))
                current_RT = current_RT(1:end-1);
            end
        end
        
        powerRTcorr_surr_metadata.f = scalogram_metadata.f;
        powerRTcorr_surr_metadata.twin = scalogram_metadata.twin;
        powerRTcorr_surr_metadata.t = scalogram_metadata.t;
        powerRTcorr_surr_metadata.eventList = scalogram_metadata.eventList;
        powerRTcorr_surr_metadata.Fs = scalogram_metadata.Fs;    
        powerRTcorr_surr_metadata.numSurrogates = numSurrogates;
        
        eventList = powerRTcorr_surr_metadata.eventList;
            
        numEvents = length(powerRTcorr_surr_metadata.eventList);
        numSamps = size(W,2);
        numFreqs = length(powerRTcorr_surr_metadata.f);
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};
            
            gabor_name = [sessionChannels{iCh}.name '_' trialType '_scalograms.mat'];
            gabor_name = fullfile(gabor_sessionDir, gabor_name);
            if ~exist(gabor_name,'file')
                disp(['file ' gabor_name ' does not exist.']);
                continue;
            end
            
            power_RTcorr_surr_name = ['power_RTcorr_surr_' sessionChannels{iCh}.name '.mat'];
            power_RTcorr_surr_name = fullfile(powerRTcorr_sessionDir, power_RTcorr_surr_name);
            
            powerRTcorr_surr = NaN(numSurrogates,numEvents, numFreqs, numSamps);
            
            % load analytic signal around each event and calculate
            % correlation coefficients
            if exist(power_RTcorr_surr_name, 'file')
                continue
            end
            
            if iCh > 1; load(gabor_name); end    % the first channel gabor scalograms were loaded earlier
            
            gaborPower = abs(W).^2;
            nTrials = length(current_RT);
            for iSurr = 1 : numSurrogates
                for iEvent = 1 : numEvents
    %                 tic
                    disp(sprintf('%s, %s, %d of %d sessions', ...
                                 ch.name, eventList{iEvent}, iSession, numSessions));

                    for iFreq = 1 : numFreqs

                        for iSamp = 1 : numSamps
                            [cc, p] = corr(squeeze(gaborPower(iEvent,iSamp,:,iFreq)), current_RT(randperm(nTrials))', ...
                                           'type','spearman');
                            powerRTcorr_surr(iSurr, iEvent, iFreq, iSamp) = cc;
                        end

                    end
    %                 toc
                end
            end    % iSurr

%             powerRTcorr_surr_mean = squeeze(nanmean(powerRTcorr_surr,1));
%             powerRTcorr_surr_std = squeeze(nanstd(powerRTcorr_surr,0,1));
%             save(power_RTcorr_surr_name, 'powerRTcorr_surr_mean', 'powerRTcorr_surr_std','powerRTcorr_surr_metadata');
            save(power_RTcorr_surr_name, 'powerRTcorr_surr','powerRTcorr_surr_metadata');
            fprintf('%s saved\n', power_RTcorr_surr_name);
            
        end
    end
    
end
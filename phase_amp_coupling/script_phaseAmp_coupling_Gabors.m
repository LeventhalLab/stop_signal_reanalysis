% script_phaseAmp_coupling_Gabors

bitOrder = 'b';
low_freq_range  = [0 100];
high_freq_range = [10 550];

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
scalogramDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/trial_scalograms';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed_Gabors';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

for i_chDB = 5 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end

    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    
    subject_scalogramDir = fullfile(scalogramDir, [implantID '_ps']);
    if ~exist(subject_scalogramDir,'dir'); continue; end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
    channels = eval( chDB_info.name );

    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );

    for iTrialType = 1 : length(trialTypeList)
        trialType = trialTypeList{iTrialType};
        
        for iCh = 1 : length(channels)
            ch = channels{iCh};
            session_scalogramDir = fullfile(subject_scalogramDir,[ch.session '_scalograms']);
            test_ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
            if ~exist(test_ch_scalogramName,'file');continue;end
            break;
        end
        if ~exist(test_ch_scalogramName,'file');continue;end
        load(test_ch_scalogramName);
        t = scalogram_metadata.t; f = scalogram_metadata.f;
        numSamples = length(t);
        phaseAmp_metadata.f = f;
        phaseAmp_metadata.t = t;
        
        numSamps  = length(t); numFreqs = length(f);
        
        phase_f_idx = find(f > low_freq_range(1) & f < low_freq_range(2));
        amp_f_idx   = find(f > high_freq_range(1) & f < high_freq_range(2));
        
        phaseAmp_metadata.phase_f   = f(phase_f_idx);
        phaseAmp_metadata.amp_f     = f(amp_f_idx);
        phaseAmp_metadata.eventList = scalogram_metadata.eventList;
        phaseAmp_metadata.trialType = trialType;
        phaseAmp_metadata.eventtWin = scalogram_metadata.twin;
        phaseAmp_metadata.Fs        = scalogram_metadata.Fs;
        
        num_fphase = length(phaseAmp_metadata.phase_f);
        num_famp   = length(phaseAmp_metadata.amp_f);
        
        numEvents = length(phaseAmp_metadata.eventList);

        if i_chDB == 7 && iTrialType == 1
            startSession = 1;
        else
            startSession = 1;
        end
        for iSession = startSession : numSessions
            session_scalogramDir = fullfile(subject_scalogramDir,[sessionList{iSession} '_scalograms']);
            if ~exist(session_scalogramDir,'file'); continue; end
            
            fprintf('trialtype: %s, session %s, %d of %d\n', trialType, sessionList{iSession}, iSession, numSessions)

            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                mkdir(phaseAmp_sessionDir);
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
            cp.locationSubClass = ROI_list;
        
            if ~isempty(strfind(trialType, 'nogo'))
                cp.task = 4;
            elseif ~isempty(strfind(trialType, 'stop'))
                cp.task = 3;
            else
                cp.task = -1;
            end
            session_chList = extractChannels( cp, channels );
            sessionChannels = channels(session_chList);

            % exclude EMG, reference channels
            cp = initChanParams();
            cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            
            cp = initChanParams();
            cp.tetrode = {'e2', 'e3', 'e02','e03'};
            sessionChannels = excludeChannels(cp, sessionChannels);
            if isempty(sessionChannels); continue; end

            numCh = length(sessionChannels);
            
            for iCh = 1 : numCh             
                
                ch = sessionChannels{iCh};
                if any(ch.wire.markedGood) == 0; continue; end
                
                disp(sprintf('%s, %d of %d sessions, %d of %d channels in session %s', ...
                    ch.name, ...
                    iSession, ...
                    numSessions, ...
                    iCh, ...
                    numCh, ...
                    sessionList{iSession}));
                
                phaseAmp_name_mrl = [ch.name '_' trialType '_phase_amp_Gabor.mat'];
                phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
                if exist(phaseAmp_name_mrl, 'file'); continue; end

                phaseAmp_metadata.chName = ch.name;
                phaseAmp_metadata.region = ch.location.name;
                
                ch_scalogramName = fullfile(session_scalogramDir,[ch.name '_' trialType '_scalograms.mat']);
                if ~exist(ch_scalogramName,'file');continue;end

                load(ch_scalogramName);

                phase_angles = angle(W(:,:,:,phase_f_idx));
                phase_angles(phase_angles < 0) = phase_angles(phase_angles < 0) + 2 * pi;
                
                sig_amp = abs(W(:,:,:,amp_f_idx));

                mrv = NaN(numEvents, num_fphase, num_famp, numSamples);
                for i_fphase = 1 : num_fphase
                    for i_famp = 1 : num_famp
                        
                        if f(phase_f_idx(i_fphase)) >= f(amp_f_idx(i_famp)); continue; end
                        
                        if numEvents > 1
                            phaseAmpfunction = squeeze(sig_amp(:,:,:,i_famp)) .* squeeze(exp(1i*phase_angles(:,:,:,i_fphase)));
                        else
                            phaseAmpfunction = (sig_amp(:,:,:,i_famp)) .* (exp(1i*phase_angles(:,:,:,i_fphase)));
                        end
                        mrv(:, i_fphase, i_famp, :) = mean(phaseAmpfunction, 3);
                        
                    end 
                end
                
%                 re_mrv = real(mrv); im_mrv = imag(mrv);
                save(phaseAmp_name_mrl, 'mrv', 'phaseAmp_metadata');

            end

        end    % for iSession...

    end   % for iTrialType...
    
end
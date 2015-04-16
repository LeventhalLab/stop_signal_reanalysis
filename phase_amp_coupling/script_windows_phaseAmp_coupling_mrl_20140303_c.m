% script_windows_phaseAmp_coupling_mrl_20140303_c
% script to calculate phase-amplitude coupling in discrete windows

bitOrder = 'b';
low_freq_range  = [0 21];
high_freq_range = [10 101];


chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventList{1} = {'noseCenterIn'};
eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
eventList{3} = eventList{2};
eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};

numTrialTypes = length(trialTypeList);
eventtWin   = zeros(numTrialTypes, 2);
% analysisWin = zeros(numTrialTypes, 1);
% stepSize    = zeros(numTrialTypes, 1);
eventtWin(1,:) = [-1 2];   % for analysis of all trials
% analysisWin(1) = 3;
% stepSize(1)    = 3;
for iTrialType = 2 : numTrialTypes
    eventtWin(iTrialType,:) = [-1 1];
%     analysisWin(iTrialType) = 0.1;
%     stepSize(iTrialType)    = 0.05;
end

for i_chDB = 7 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end

    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);
    subject_hilbertDir_025Hz = fullfile(hilbert_025Hz_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir_1Hz, 'dir')
        disp([subject_hilbertDir_1Hz ' not found. Skipping ' implantID '...'])
        continue
    end
    if ~exist(subject_hilbertDir_025Hz, 'dir')
        disp([subject_hilbertDir_025Hz ' not found. Skipping ' implantID '...'])
        continue
    end

    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end

    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );

    if i_chDB == 7
        startTrialType = 1;
    else
        startTrialType = 1;
    end
    for iTrialType = startTrialType : length(trialTypeList)
        trialType = trialTypeList{iTrialType}
        
        phaseAmp_metadata.low_freq_range  = low_freq_range;
        phaseAmp_metadata.high_freq_range = high_freq_range;
        phaseAmp_metadata.eventList       = eventList{iTrialType};
        phaseAmp_metadata.trialType       = trialType;
        phaseAmp_metadata.eventtWin       = eventtWin(iTrialType, :);
%         phaseAmp_metadata.analysisWin     = analysisWin(iTrialType);
%         phaseAmp_metadata.stepSize        = stepSize(iTrialType);

        if i_chDB == 7 && iTrialType == 1
            startSession = 1;
        else
            startSession = 1;
        end
        for iSession = startSession : numSessions
            
            if strcmpi(sessionList{iSession},'IM328_20120914_08-10-23');continue;end
            
            disp(sprintf('trialtype: %s, session %s, %d of %d', trialType, sessionList{iSession}, iSession, numSessions))

            phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession}, [sessionList{iSession} '_' trialType]);
            if ~exist(phaseAmp_sessionDir, 'dir')
                mkdir(phaseAmp_sessionDir);
            end
            
            cp = initChanParams();
            cp.session = sessionList{iSession};
        
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

            hilbert_sessionDir_1Hz = fullfile(subject_hilbertDir_1Hz, sessionList{iSession});
            hilbert_sessionDir_025Hz = fullfile(subject_hilbertDir_025Hz, sessionList{iSession});

            metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
            metadata_filename_1Hz = fullfile(hilbert_sessionDir_1Hz, metadata_filename);
            metadata_filename_025Hz = fullfile(hilbert_sessionDir_025Hz, metadata_filename);

            if ~exist(metadata_filename_1Hz, 'file')
                continue
    %             error([metadata_filename ' could not be found.']);
            end
            if ~exist(metadata_filename_025Hz, 'file')
                continue
    %             error([metadata_filename ' could not be found.']);
            end
            md_1Hz   = load(metadata_filename_1Hz);
            md_025Hz = load(metadata_filename_025Hz);
            
            centerFreqs_1Hz   = mean(md_1Hz.metadata.freqBands, 2);
            centerFreqs_025Hz = mean(md_025Hz.metadata.freqBands, 2);
            
            max025Hz = max(centerFreqs_025Hz);

            low_freq_idx_025  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
            low_freq_idx_1    = find(centerFreqs_1Hz >= max025Hz & centerFreqs_1Hz <= low_freq_range(2));
            low_freq_idx      = [low_freq_idx_025; low_freq_idx_1];
            high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

            low_freqs  = [centerFreqs_025Hz(low_freq_idx_025); centerFreqs_1Hz(low_freq_idx_1)];
            high_freqs = centerFreqs_1Hz(high_freq_idx);
            
            num_low_freq_025 = length(low_freq_idx_025);
            num_low_freq  = length(low_freq_idx);
            num_high_freq = length(high_freq_idx);

            phaseAmp_metadata.Fs = md_025Hz.metadata.Fs;

            for iCh = 1 : numCh             
                
                ch = sessionChannels{iCh};
                if any(ch.wire.markedGood) == 0; continue; end
                
                if strcmp(ch.name, 'IM296_20120402_10-57-37t14') || ...
                   strcmp(ch.name, 'IM296_20120405_11-50-15t04') || ...
                   strcmp(ch.name, 'IM296_20120405_11-50-15t05') || ...
                   strcmp(ch.name, 'IM296_20120406_10-17-28t04') || ...
                   strcmp(ch.name, 'IM296_20120406_10-17-28t05') || ...
                   strcmp(ch.name, 'IM296_20120408_14-00-40t04') || ...
                   strcmp(ch.name, 'IM296_20120408_14-00-40t05') || ...
                   strcmp(ch.name, 'IM296_20120409_10-46-10t05') || ...
                   strcmp(ch.name, 'IM296_20120410_10-01-04t02') || ...
                   strcmp(ch.name, 'IM296_20120410_10-01-04t04') || ...
                   strcmp(ch.name, 'IM296_20120410_10-01-04t05') || ...
                   strcmp(ch.name, 'IM296_20120411_10-15-38t04') || ...
                   strcmp(ch.name, 'IM296_20120411_10-15-38t05') || ...
                   strcmp(ch.name, 'IM296_20120412_09-58-01t04') || ...
                   strcmp(ch.name, 'IM296_20120412_09-58-01t05') || ...
                   strcmp(ch.session, 'IM310_20120516_10-18-04')
                    continue
                end
                
                disp(sprintf('%s, %d of %d sessions, %d of %d channels in session %s', ...
                    ch.name, ...
                    iSession, ...
                    numSessions, ...
                    iCh, ...
                    numCh, ...
                    sessionList{iSession}));
                
                phaseAmp_name_mrl = [ch.name '_' trialType '_phase_amp_mrl.mat'];
                phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
                if exist(phaseAmp_name_mrl, 'file'); continue; end

                phaseAmp_metadata.chName = ch.name;
                phaseAmp_metadata.region = ch.location.name;

                [mrv, low_freqs, high_freqs] = phase_amp_Canolty_mrl_20140324(ch, ...
                                                                              eventList{iTrialType}, ...
                                                                              eventtWin(iTrialType, :), ...
                                                                              trialType, ...
                                                                              'lowfreqrange', low_freq_range, ...
                                                                              'highfreqrange', high_freq_range, ...
                                                                              'hilbert1Hzdir', hilbert_1Hz_directory, ...
                                                                              'hilbert025Hzdir', hilbert_025Hz_directory );
                % the function mean_phase_amp_comodugram calculates
                % phase-amplitude coupling based on the algorithm of Tort et
                % al, PNAS 2008. It looks at corresponding windows around
                % different behavioral events, and averages the phase-amplitude
                % comodugram across trials.
                
                re_mrv = real(mrv); im_mrv = imag(mrv);
                phaseAmp_metadata.low_freqs = low_freqs; phaseAmp_metadata.high_freqs = high_freqs;
                save(phaseAmp_name_mrl, 're_mrv', 'im_mrv', 'phaseAmp_metadata');

            end

        end    % for iSession...

    end   % for iTrialType...
    
end
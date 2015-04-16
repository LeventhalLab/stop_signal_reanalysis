% script_phaseAmp_surrogates_c3

bitOrder = 'b';
low_freq_range  = [0 10];
high_freq_range = [10 101];


chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_1Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
hilbert_025Hz_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins';
phaseRThist_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_histogram_analysis';
phaseAmp_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phaseAmp_windowed';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

eventLists{1} = {'noseCenterIn'};
% eventList{2} = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
% eventList{3} = eventList{2};
% eventList{4} = {'cueOn','noseCenterIn','tone','whiteNoise','foodHopperClick'};
% eventList{5} = {'cueOn','noseCenterIn','tone','whiteNoise','noseCenterOut'};
% eventList{6} = {'cueOn','noseCenterIn','whiteNoise','foodHopperClick'};
% eventList{7} = {'cueOn','noseCenterIn','whiteNoise','noseCenterOut'};
% trialTypeList = {'any','correctgo', 'wronggo', 'correctstop', 'failedstop', 'correctnogo', 'failednogo'};
trialTypeList = {'any'};
trialType = trialTypeList{1};
iTrialType = 1;

numTrialTypes = length(trialTypeList);
eventtWin   = [-1,2];
analysisWin = 3;
stepSize    = 3;

numSurrogates = 200;
maxSkip = 100;    % in seconds
minSkip = 10;

for i_chDB = 3 : 3%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end

    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
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
    subject_hilbertDir_1Hz = fullfile(hilbert_1Hz_directory, [implantID '_hilbert']);

    subject_phaseAmpdir = fullfile(phaseAmp_directory, [implantID '_phase_amp']);
    if ~exist(subject_phaseAmpdir, 'dir')
        mkdir(subject_phaseAmpdir);
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );


        
    surrogate_phaseAmp_metadata.low_freq_range  = low_freq_range;
    surrogate_phaseAmp_metadata.high_freq_range = high_freq_range;
    surrogate_phaseAmp_metadata.eventList       = eventLists{1};
    surrogate_phaseAmp_metadata.trialType       = trialType;
    surrogate_phaseAmp_metadata.eventtWin       = eventtWin(iTrialType, :);
    surrogate_phaseAmp_metadata.analysisWin     = analysisWin(iTrialType);
    surrogate_phaseAmp_metadata.stepSize        = stepSize(iTrialType);

    for iSession = 4 : numSessions

        phaseAmp_sessionDir = fullfile(subject_phaseAmpdir, sessionList{iSession});
        if ~exist(phaseAmp_sessionDir, 'dir')
            mkdir(phaseAmp_sessionDir);
        end
        
        cp = initChanParams();
        cp.session = sessionList{iSession};

        session_chList = extractChannels( cp, channels );
        sessionChannels = channels(session_chList);

        % exclude EMG, reference channels
        cp = initChanParams();
        cp.locationSubClass = {'EMG', 'EEGLAM', 'Ref'};
        sessionChannels = excludeChannels(cp, sessionChannels);
            
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

        low_freq_idx  = find(centerFreqs_025Hz >= low_freq_range(1) & centerFreqs_025Hz <= low_freq_range(2));
        high_freq_idx = find(centerFreqs_1Hz >= high_freq_range(1) & centerFreqs_1Hz <= high_freq_range(2));

        low_freq = centerFreqs_025Hz(low_freq_idx);
        high_freq = centerFreqs_1Hz(high_freq_idx);

        surrogate_phaseAmp_metadata.low_freq = low_freq;
        surrogate_phaseAmp_metadata.high_freq = high_freq;
            
        num_low_freq  = length(low_freq_idx);
        num_high_freq = length(high_freq_idx);

        surrogate_phaseAmp_metadata.Fs = md_025Hz.metadata.Fs;
%             surrogate_phaseAmp_metadata.chList = cell(1, numCh);
%             surrogate_phaseAmp_metadata.regionList = cell(1, numCh);
            
        totSamps = md_1Hz.metadata.numSamps;
            
        minSkipSamps = round(minSkip * surrogate_phaseAmp_metadata.Fs);
        maxSkipSamps = round(maxSkip * surrogate_phaseAmp_metadata.Fs);
        skipSampRange = maxSkipSamps - minSkipSamps;
        skip = minSkipSamps + ceil(skipSampRange .* rand(numSurrogates,1));
%             skip(find(skip>maxSkip))=[];
%             skip(find(skip<minSkip))=[];
            
        for iCh = 1 : numCh
            iCh
            ch = sessionChannels{iCh};
                
            phaseAmp_name_mrl = [ch.name '_surr_fullSig_phase_amp_mrl.mat'];
            phaseAmp_name_mrl = fullfile(phaseAmp_sessionDir, phaseAmp_name_mrl);
            if exist(phaseAmp_name_mrl, 'file')
                load(phaseAmp_name_mrl);
                if surrogate_phaseAmp_metadata.numLowFreqsComplete == num_low_freq
                    continue
                end
                numLowFreqsComplete = surrogate_phaseAmp_metadata.numLowFreqsComplete;
            else
                numLowFreqsComplete = 0;
            end
                
            surrogate_phaseAmp_metadata.chName = ch.name;
            surrogate_phaseAmp_metadata.region = ch.location.name;

            as2 = zeros(num_high_freq, totSamps);
            as1 = zeros(num_low_freq, totSamps);
            
            hilbert_name_1Hz = ['analytic_' ch.name '.bin'];
            hilbert_name_1Hz = fullfile(hilbert_sessionDir_1Hz, hilbert_name_1Hz);
            
            hilbert_name_025Hz = ['analytic_' ch.name '.bin'];
            hilbert_name_025Hz = fullfile(hilbert_sessionDir_025Hz, hilbert_name_025Hz);

            for i_f2 = 1 : num_high_freq
                as2(i_f2,:) = readAnalyticSignal(hilbert_name_1Hz, md_1Hz.metadata, [0, md_1Hz.metadata.duration], i_f2);
            end

%                         temp = getAnalyticAroundEvent_20140221( ch, ...
%                                                                 high_freq_idx(i_f2), ...
%                                                                 eventList{iEventType}, ...
%                                                                 trialType, ...
%                                                                 eventtWin, ...
%                                                                 'hilbertdir', hilbert_1Hz_directory );
%                         as2(i_f2, :) = reshape(temp', 1, totSamps);
%                     end
                    
%                     as2 = squeeze(reshape(as2, num_high_freq, numTrials * sampsPerTrial, 1));
            sig_amp = abs(as2);

            for i_f1 = numLowFreqsComplete + 1  : num_low_freq
                as1(i_f1,:) = readAnalyticSignal(hilbert_name_025Hz, md_025Hz.metadata, [0, md_025Hz.metadata.duration], i_f1);
            end

%                         temp = getAnalyticAroundEvent_20140221( ch, ...
%                                                                 low_freq_idx(i_f1), ...
%                                                                 eventList{iEventType}, ...
%                                                                 trialType, ...
%                                                                 eventtWin, ...
%                                                                 'hilbertdir', hilbert_025Hz_directory );
%                                                             
%                         as1(i_f1, :) = reshape(temp', 1, totSamps);
% 
%                     end
                    
%                     as1 = squeeze(reshape(as1, num_high_freq, numTrials * sampsPerTrial, 1));
            phase_angles = angle(as1);
            surrogate_mean = zeros(num_low_freq, num_high_freq);
            surrogate_std  = zeros(num_low_freq, num_high_freq);
            for i_f1 = numLowFreqsComplete + 1 : num_low_freq
                i_f1
                for i_f2 = 1 : num_high_freq

                    if centerFreqs_025Hz(low_freq_idx(i_f1)) >= centerFreqs_1Hz(high_freq_idx(i_f2)); continue; end

                    amp_f2 = sig_amp(i_f2, :);
                    surrogate_mrv = zeros(1, numSurrogates);
                    for iSurrogate = 1 : numSurrogates
%                                 surrogate_amp(iSurrogate, :) = [squeeze(sig_amp(i_f2, :, skip(iSurrogate):end)), ...
%                                                                    squeeze(sig_amp(i_f2, :, 1:skip(iSurrogate)-1))];
                        surrogate_amp = circshift(amp_f2, [0, skip(iSurrogate)]);
                        surrogate_mrv(iSurrogate) = abs(mean(surrogate_amp.*exp(1i*squeeze(phase_angles(i_f1, :)))));
                    end    % for iSurrogate...
                    surrogate_mean(i_f1,i_f2) = mean(surrogate_mrv);
                    surrogate_std(i_f1, i_f2) = std(surrogate_mrv);
%                             [surrogate_mean(i_f1, i_f2), surrogate_std(i_f1, i_f2)] = normfit(surrogate_mrl);

                end    % for i_f2...
                        
                numLowFreqsComplete = numLowFreqsComplete + 1;
                surrogate_phaseAmp_metadata.numLowFreqsComplete = numLowFreqsComplete;
                save(phaseAmp_name_mrl, 'surrogate_mean', 'surrogate_std', 'surrogate_phaseAmp_metadata');

            end    % for i_f1...
                
        end    % for iCh = 1 : numCh
    end    % for iSession...
        
end    % for i_chDB...
                
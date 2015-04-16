% script_storeHilbertPowerComod

bitOrder = 'b';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
comod_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power-power comodugrams';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

for i_chDB = 1 : 2: length(chDB_list)
    
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
    
    subject_comoddir = fullfile(comod_directory, [implantID '_phase-power_comods']);
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
        cp.locationSubClass = {'EMG', 'EEGLAM'};
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
            error([metadata_filename ' could not be found.']);
        end
        load(metadata_filename);
        freqList = mean(metadata.freqBands, 2);
        
        for iCh1 = 1 : numCh
            ch1 = sessionChannels{iCh1};
            tetName1 = ch1.name(end-2:end);
            
            hilbert_name1 = ['analytic_' sessionChannels{iCh1}.name '.bin'];
            hilbert_name1 = fullfile(hilbert_sessionDir, hilbert_name1);

            
            for iCh2 = iCh1 : numCh
                ch2 = sessionChannels{iCh2};
                tetName2 = ch2.name(end-2:end);
                
                hilbert_name2 = ['analytic_' sessionChannels{iCh2}.name '.bin'];
                hilbert_name2 = fullfile(hilbert_sessionDir, hilbert_name2);
            
                comod_name = [sessionList{iSession} '_' tetName1 '_' tetName2 '_power_power_comodugram.mat'];
                comod_name = fullfile(comod_sessionDir, comod_name);
            
                if exist(comod_name, 'file')
                    continue;
                end
                
                ppcomod = zeros(length(freqList));
        
                for if1 = 1 : length(freqList)
                    if1
                    as1 = readAnalyticSignal(hilbert_name1, metadata, ...
                                             [0 metadata.duration], ...
                                             if1);
                    for if2 = 1 : length(freqList)
                        as2 = readAnalyticSignal(hilbert_name2, metadata, ...
                                                 [0 metadata.duration], ...
                                                 if2);
                        temp = corrcoef(abs(as1).^2, abs(as2).^2);
                        ppcomod(if1, if2) = temp(1,2);
                        
                    end
                end
                                                           
                save(comod_name, 'ppcomod', 'freqList', 'metadata');
                
            end
            
        end
        
    end
    
end
            
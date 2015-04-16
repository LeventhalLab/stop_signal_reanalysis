% script to call calc_SFC_freqbands

[chDB_list, chDB_fnames, spikeDB_list, spikeDB_fnames] = get_chStructs_for_analysis();

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
saveRoot = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Fries_sfc';

Fsfc_metadata.eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn'};
Fsfc_metadata.trialType = 'correctgo';
Fsfc_metadata.eventWin = [-1 1];
Fsfc_metadata.spikeWin = 0.1;
Fsfc_metadata.stepSize = 0.02;   % 100 ms windows in which to look for spikes
Fsfc_metadata.sta_win  = 3;      % multiples of a single period at a given frequency within which to calculate the sta.
                   % that is, this is the duration of each spike-triggered
                   % lfp signal in number of full-wave periods
numSteps = floor(range(Fsfc_metadata.eventWin) / Fsfc_metadata.stepSize);
numEventTypes = length(Fsfc_metadata.eventList);
numFreqs = 99;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant spike and channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    if ~exist(spikeDB_list{i_chDB}, 'var')
        spikeDB_file = fullfile(chDB_directory, spikeDB_fnames{i_chDB});
        disp(['loading ' spikeDB_file]);
        load( spikeDB_file );
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
        spikeDB_info = whos( [spikeDB_list{i_chDB}(1:3) 'spike'] );
    else
        implantID = chDB_list{i_chDB}(1:5);
        chDB_info = whos( [implantID 'Ch*'] );
        spikeDB_info = whos( implantID );
    end
    channels = eval( chDB_info.name );
    spike    = eval( spikeDB_info.name );
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    
    subject_saveDir = fullfile(saveRoot, [implantID, '_Fries_sfc']);
    if ~exist(subject_saveDir, 'dir')
        mkdir(subject_saveDir);
    end
    
    sessionList = getSessionsfromChannelDB( spike );
    
    for iSession = 1 : length(sessionList)
        
        cp = initChanParams();
        cp.session = sessionList{iSession};
        disp(sessionList{iSession});
        
        saveName = [implantID, '_', sessionList{iSession}, '.mat'];
        saveName = fullfile(subject_saveDir, saveName);
        if exist(saveName,'file'); continue; end
        
        spikeList = extractChannels(cp, spike);
        sessionSpikes = spike(spikeList);
        
        numSessionUnits = length(sessionSpikes);
        sfc = zeros(numSessionUnits, numEventTypes, numSteps, numFreqs);
        numSpikes = zeros(numSessionUnits, numEventTypes, numSteps);
        
        spName = cell(1, numSessionUnits);
        for iUnit = 1 : numSessionUnits
            disp(sprintf('unit %d of %d', iUnit, numSessionUnits));
            
            sp = sessionSpikes{iUnit};
            spName{iUnit} = sp.name;
            
            cp = initChanParams();
            cp.tetrode = sp.tetrode.name;
            cp.session = sp.session;
            
            chList = extractChannels(cp, channels);
            if isempty(chList); continue; end
            
            ch = channels{chList};
            sp_ts = sp.timestamps.spike;
            
            [unit_sfc, unit_numSpikes, f, numTrials] = ...
                calc_SFC_freqbands( ch, sp_ts, Fsfc_metadata.trialType, ...
                Fsfc_metadata.eventList, Fsfc_metadata.eventWin, Fsfc_metadata.spikeWin, ...
                Fsfc_metadata.stepSize, ...
               'hilbert_directory', hilbert_directory, 'stawin', sta_win );
           
            sfc(iUnit, :, :, :) = unit_sfc;
            numSpikes(iUnit, :, :) = unit_numSpikes;
                           
        end
        Fsfc_metadata.f = f;
        metadata.numTrials = numTrials;
        
        save(saveName, 'sfc', 'numSpikes', 'Fsfc_metadata','spName');
        
    end
            
end
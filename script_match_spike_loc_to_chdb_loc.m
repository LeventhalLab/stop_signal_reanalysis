% script to copy location labels from Nico and Robert's spike DBs into the
% channel DB

% WHY DON'T NUMBER OF TRIALS MATCH BETWEEN THE SPIKE AND CHANNEL DATABASES?

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
[chDB_list, chDB_fnames, spikeDB_list, spikeDB_fnames] = ...
    get_chStructs_for_analysis;

for i_chDB = 5 : length(chDB_list)
    
    % first, load the relevant channel and spike DBs, if necessary
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
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
        spikeDB_info = whos( [spikeDB_list{i_chDB}(1:3) 'spike'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
        spikeDB_info = whos( implantID );
    end
    channels = eval( chDB_info.name );
    spike = eval( spikeDB_info.name );
    
    num_spikeCh = length(spike);
    % match tetrodes based on date, tetrode number, and number of trials
    % (should exclude sleep sessions)
    for i_spike = 1 : num_spikeCh
        sp = spike{i_spike};
        
        cp = initChanParams();
        
        % extract the date from the session name in the spike DB
        spikeDate = sp.session(7:14);
        spikeDateVector = datevec(spikeDate,'yyyymmdd');
        chDateStr = datestr(spikeDateVector, 'yyyy-mm-dd');
        cp.date = chDateStr;
        cp.tetrode = sp.tetrode.name(1:3);
        
        chList = extractChannels(cp, channels);
        if isempty(chList); continue; end
        matchChannel = [];
        if length(chList) > 1
            sp_hsdName = sp.hsdFilename(1:end-4);
            for iCh = 1 : length(chList)
                if strcmpi(channels{chList(iCh)}.session,sp_hsdName)
                    matchChannel = chList(iCh);
                    break;
                end
            end
        else
            matchChannel = chList;
        end
        if isempty(matchChannel); continue; end
        
        channels{matchChannel}.location.name = sp.location.name;
        channels{matchChannel}.location.subclass = sp.location.subClass;
        channels{matchChannel}.location.ml = sp.location.ml;
        channels{matchChannel}.location.ap = sp.location.ap;
        channels{matchChannel}.location.dv = sp.location.dv;
    end    % for i_spike = 1 : num_spikeCh
    
    eval( [chDB_info.name ' = channels;'] );
    
    new_chDBname = [chDB_fnames{i_chDB}(1:end-4) '_loc.mat'];
    save(fullfile(chDB_directory, new_chDBname), chDB_info.name);
    
end    % for i_chDB ...
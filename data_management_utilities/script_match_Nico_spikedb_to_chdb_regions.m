chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

for i_chDB = 5 : 5%length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:5));
    
    % load the corresponding spike db, if necessary
    if ~exist(implantID, 'var')
        chDB_file = fullfile(chDB_directory, [implantID '.mat']);
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    spikedb = eval(implantID);
    
    chDB_info = whos( [chDB_list{i_chDB}(1:5) 'Ch*'] );
    chdb      = eval( chDB_info.name );
    
    for iSpike = 1 : length(spikedb)
        
        
    
    
    
    
end
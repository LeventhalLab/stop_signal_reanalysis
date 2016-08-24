% script_extractLFPs

targetFs = 1000;
rootReadDir = '/Volumes/Leventhal_lab_HD01';
% rootReadDir = '/Volumes/Untitled';
rootSaveDir = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/high_cutoff_stop-signal LFPs';
chDB_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';

[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);

for i_chDB = 11 : 11%length(chDB_list)

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
    
    subject_HSDDir = fullfile(rootReadDir, [implantID '-rawdata']);
    if ~exist(subject_HSDDir,'dir'); continue; end
    
    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
    channels = eval( chDB_info.name );
    
    subject_LFPDir = fullfile(rootSaveDir, [implantID '_HF_LFPs']);
    if ~exist(subject_LFPDir, 'dir')
        mkdir(subject_LFPDir);
    end
    
    sessionList = getSessionsfromChannelDB( channels );
    numSessions = length( sessionList );
    
%     for iSession = 1 : numSessions
        
%         cp = initChanParams();
%         cp.session = sessionList{iSession};
%         cp.locationSubClass = ROI_list;
%         
%         sessionChannels = channels(extractChannels( cp, channels ));
%         
%         if isempty(sessionChannels); continue; end
%         sessionDate = sessionList{iSession}(4:end);
%         
%         session_HSDDir = fullfile(subject_HSDDir,[implantID '_' sessionDate 'a']);
%         if ~exist(session_HSDDir,'dir'); continue; end
        
%         session_LFPName = fullfile(subject_LFPDir, [implantID '_' sessionDate);
%         if ~exist(phaseAmp_sessionDir, 'dir')
%             mkdir(phaseAmp_sessionDir);
%         end

        if i_chDB < 5
            extractLFPbyChannels_20151204(channels, ...
                                         'targetfs', targetFs, ...
                                         'savedir', subject_LFPDir, ...
                                         'opendir', subject_HSDDir);
        else
            extractLFPbyChannels_20151204_Nico(channels, ...
                                         'targetfs', targetFs, ...
                                         'savedir', subject_LFPDir, ...
                                         'opendir', subject_HSDDir);
        end
        
%     end
    
end
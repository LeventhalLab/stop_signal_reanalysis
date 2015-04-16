% script_create_phaseRTcorr_summaries

% hilbert_directory     = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
phaseRTcorr_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

for i_chDB = 1 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    
    subject_phaseRTcorrdir = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr']);
    if ~exist(subject_phaseRTcorrdir, 'dir')
        continue;
    end
    
    chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    channels = eval( chDB_info.name );
    
    sessionList = getSessionsfromChannelDB(channels);
    numSessions = length(sessionList);
    
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
        
        phaseRTcorr_sessionDir = fullfile(subject_phaseRTcorrdir, sessionList{iSession});
        if ~exist(phaseRTcorr_sessionDir, 'dir')
            continue;
        end
        
        for iCh = 1 : numCh
            
            ch = sessionChannels{iCh};

            surrogateName = ['phase_RT_surrogates_' ch.name '.mat'];
            surrogateName = fullfile(phaseRTcorr_sessionDir, surrogateName);
            if ~exist(surrogateName, 'file'); continue; end
            
            surrogateMetadataName = ['phase_RT_surrogate_metadata_' ch.name '.mat'];
            surrogateMetadataName = fullfile(phaseRTcorr_sessionDir, surrogateMetadataName);
            if ~exist(surrogateMetadataName, 'file'); continue; end
            
            phase_RTcorr_name = ['phase_RT_analysis_' sessionChannels{iCh}.name '.mat'];
            phase_RTcorr_name = fullfile(phaseRTcorr_sessionDir, phase_RTcorr_name);
            if ~exist(phase_RTcorr_name, 'file'); continue; end
            
            load(surrogateName);
            load(surrogateMetadataName);
            load(phase_RTcorr_name);
            
            h_fig = phaseRTcorr_summary_sheet(ch.name, phaseRTcorr_metadata, RTphases, surrogate_metadata, mrl);
            
            summaryName = ['phaseRTcorr_plots_', ch.name, '.pdf'];
            summaryName = fullfile(phaseRTcorr_sessionDir, summaryName);
            
            for iFig = 1 : length(h_fig)
                
                if iFig == 1
                    figure(h_fig(iFig))
                    export_fig(summaryName,'-pdf','-q101','-painters')
                else
                    figure(h_fig(iFig))
                    export_fig(summaryName,'-pdf','-q101','-painters','-append')
                end
            end
        end
    end
end
                    
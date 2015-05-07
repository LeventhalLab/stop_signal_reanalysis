% script to take calculated power spectra, find the important specral peaks
% by whitening the spectra, then find the width of the peaks

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
powerSpectraDir = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/power_spectra';
lfp_root          = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs';
sample_lfp_name = '/Volumes/PublicLeventhal2/dan/stop-signal reanalysis/stop-signal LFPs/IM164_LFPs/IM164_20091112_13-56-32.hsdf';

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

ROI_list = {'eegorb','cpu','gp','stn','snr'};
numRegions = length(ROI_list);
flength = 1025;    % just to preallocate matrices


ylim = [-15 0];
max_f_to_plot = [10,100];
for i_chDB = 1:4%length(chDB_list)
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
    else
        implantID = chDB_list{i_chDB}(1:5);
    end
    
    subject_powerSpectraDir = fullfile(powerSpectraDir, [implantID '_ps']);
    if ~exist(subject_powerSpectraDir,'dir')
        continue
    end
    PDFname = fullfile(subject_powerSpectraDir, [implantID '_ps_plots.pdf']);
    
    cd(subject_powerSpectraDir);
    fileList = dir('ps_*');
    numSessions = length(fileList);
    
    
    mean_session_psd = NaN(numSessions, numRegions, flength);
    numPages = 0;
    for iSession = 1 : numSessions

        fname = fileList(iSession).name;
        load(fname);
        
        numCh = length(ps_metadata.region);
        for iRegion = 1 : numRegions
            
            numRegionChannels = 0;regionChList = [];
            for iCh = 1 : numCh
                if strcmpi(ps_metadata.region{iCh},ROI_list{iRegion})
                    numRegionChannels = numRegionChannels + 1;
                    regionChList = [regionChList,iCh];
                end
            end
            if numRegionChannels > 0
                session_ch_means = NaN(numRegionChannels, flength);

                for iRegionCh = 1 : numRegionChannels
                    session_ch_means(iRegionCh,:) = mean(squeeze(pxx(regionChList(iRegionCh),:,:)),2)';
                end
                mean_session_psd(iSession, iRegion, :) = nanmean(session_ch_means,1);
            end
            
        end   % for iRegion...
            
    end    % for iSession...
    

    % mean across sessions:    plot(f(f<max_f_to_plot(ii)), log(squeeze(nanmean(mean_session_psd(:,iRegion,f<max_f_to_plot(ii)), 1))),'r')
    for iRegion = 1 : numRegions
        test_signal = smooth(squeeze(nanmean(mean_session_psd(:,iRegion,:),1)).*f,5);
        [peakLocations, peakValues] = peakdetect(test_signal);
    end
    
%     spectral_peaks{i_chDB,i_Region} = 
end    % for i_chDB...
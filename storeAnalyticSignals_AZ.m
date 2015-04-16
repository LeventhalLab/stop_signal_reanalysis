% script to calculate and store the analytic signal for all relevant field
% potentials

bitOrder = 'b';

chDB_directory    = '\\141.214.46.83\PublicLeventhal1\dan\stop-signal reanalysis\stop-signal data structures';
lfp_root          = '\\141.214.46.83\PublicLeventhal1\dan\stop-signal reanalysis\stop-signal LFPs';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

% establish the center frequencies for the band pass filters. Make the
% lowest low-cut boundary 1 Hz.
freq_band_width = 1;
center_freq     = ((1 + (freq_band_width/2)) : freq_band_width : 100)';
numBands        = length(center_freq);
band_edge_width = 0.5;
% freq_bands = 1 : 2 : 101;    % borders of frequency bands at which to calculate the hilbert transforms

hilbert_directory = sprintf('\\141.214.46.83\PublicLeventhal1\dan\stop-signal reanalysis\Hilbert transformed LFP %d Hz bins', freq_band_width);
hilbert_directory = '\\141.214.46.83\PublicLeventhal1\dan\stop-signal reanalysis\Hilbert transformed LFP 1 Hz bins';

% get frequency bands
M = [0 1 0];
freqBands = zeros(numBands, 4);
freqBands(:,1) = (center_freq - freq_band_width/2) - band_edge_width;
freqBands(:,2) = (center_freq - freq_band_width/2);
freqBands(:,3) = (center_freq + freq_band_width/2);
freqBands(:,4) = (center_freq + freq_band_width/2) + band_edge_width;
    
for i_chDB = 3 : 3%length(chDB_list)
    
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
    subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);
    if ~exist(subject_hilbertDir, 'dir')
        mkdir(subject_hilbertDir);
    end
    
    lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
    cd(lfp_directory);

    if i_chDB < 5
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        chDB_info = whos( [implantID 'Ch*'] );
    end
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
        cp.locationSubClass = {'EMG', 'EEGLAM', 'REF'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        cp = initChanParams();
        cp.tetrode = {'e2', 'e3'};
        sessionChannels = excludeChannels(cp, sessionChannels);
        if isempty(sessionChannels); continue; end
        
        numCh = length(sessionChannels);
        
        % load the header for the current session
        sessionDate = sessionChannels{1}.date;
        if length(sessionDate) > 8
            sessionDate = datestr(sessionDate, 'yyyymmdd');
        end
        
        lfp_fileinfo = dir(['*_' sessionDate '*.hsdf']);
        if isempty(lfp_fileinfo); continue; end
        
%         if length(lfp_fileinfo) ~= 1
%             disp([num2str(length(lfp_fileinfo)) ' files found for sessions on ' sessionDate]);
%             break;
%         end
        
        spikeInterp_idx = []; orig_lfp_idx = [];
        if length(lfp_fileinfo) > 1
            for iFile = 1 : length(lfp_fileinfo)
                if isempty(strfind(lfp_fileinfo(iFile).name, 'spikeInterp'))    % use only original files, ignore the "spike interpolation" files
                    orig_lfp_idx = [orig_lfp_idx, iFile];
                else
                    spikeInterp_idx = [spikeInterp_idx, iFile];
                end
            end
            lfp_fileinfo = lfp_fileinfo(orig_lfp_idx);
            if length(lfp_fileinfo) > 1
                disp([num2str(length(lfp_fileinfo)) ' files found for sessions on ' sessionDate]);
                break;
            end
        end
        
        lfp_fileName = lfp_fileinfo.name;
        
        % read in the LFP header and data
        header = getHSDHeader( lfp_fileName );
        numRecordingSites = header.main.num_channels;
        lfp_wireNums = zeros(1, numRecordingSites);
        for i_site = 1 : numRecordingSites
            lfp_wireNums(i_site) = header.channel(i_site).original_number;
        end
        
        Fs = lfpFs( header );
        lfpDuration = getHSDlength( 'filename', lfp_fileName );
        
        hilbert_sessionDir = fullfile(subject_hilbertDir, sessionList{iSession});
        if ~exist(hilbert_sessionDir, 'dir')
            mkdir(hilbert_sessionDir);
        end
        
        % create a text file with metadata for the files that will be written
        % to this folder
        metadata_filename = [sessionList{iSession} 'hilbert_metadata.mat'];
        metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);
        
        if exist(metadata_filename, 'file')
            load(metadata_filename);
        else
            metadata.session   = sessionList{iSession};
            metadata.duration  = lfpDuration;
            metadata.numSamps  = 0;
            metadata.sigmax    = zeros(1, numCh); 
            metadata.Fs        = Fs;
            metadata.chNames   = cell(1, numCh);
            metadata.freqBands = freqBands;
            metadata.bitOrder  = bitOrder;
            for iCh = 1 : numCh
                metadata.chNames{iCh} = sessionChannels{iCh}.name;
            end

            metadata.numWrittenChannels = 0;   % number of channels for which all frequencies have been filtered, hilbert transformed, and stored
            metadata.numWrittenFreqs    = 0;

            metadata.b = cell(1, numBands);
            metadata.comment = ['Analytic signals calculated by bandpass filtering LFPs and calculating the Hilbert transform. ' ...
                                'Data are stored in binary files as double precision complex numbers. Each analytic signal is ' ...
                                'numSamps values long. Analytic signals are stored column-wise, with the left-most column representing ' ...
                                'the LFP filtered in the freqBands(:,1) band, etc. The complex numbers are stored as 2 sets of double ' ...
                                'precision numbers - first, the full stream of real parts, then the full stream of imaginary parts '...
                                '(that is, [1+2i, 3+4i] would be stored as [1,3,2,4]. Double precision numbers were scaled to the full ' ...
                                'int16 range and written as int16s. That is if x is a 16-bit integer read from the file, to get the '...
                                'correct value, take x * sigmax(channel_index) / 32767. sigmax is an array stored in the metadata ' ...
                                'structure.'];
            save(metadata_filename, 'metadata');
        end
        
        % check to see if this session is already done
        if (metadata.numWrittenFreqs == numBands) && (metadata.numWrittenChannels == numCh)
            continue;
        end
        
        disp(['loading ' lfp_fileName '...']);
        tic
        lfp = readHSD( lfp_fileName, ...
                       numRecordingSites, ...
                       header.dataOffset, ...
                       Fs, ...
                       [0, lfpDuration] );
        toc

        numSamples = size(lfp, 2);
        metadata.numSamps = numSamples;

        t = linspace(1/Fs, lfpDuration, numSamples);
        
        if metadata.numWrittenChannels == numCh; continue; end   % if all channels for this session already written, move to next session (not sure if this is really necessary)
        chtic = tic;
        for iCh = metadata.numWrittenChannels + 1 : numCh
            iCh
            
            repWire = getRepWire( sessionChannels{iCh} );
            
            lfp_idx = find( lfp_wireNums == repWire );
            
            volt_lfp = int2volt(lfp(lfp_idx, :), 'gain', header.channel(lfp_idx).gain);

            hilbert_name = ['analytic_' sessionChannels{iCh}.name '.bin'];
            hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);
            
            if exist(hilbert_name,'file')
                fid = fopen(hilbert_name, 'a', bitOrder);
            else
                fid = fopen(hilbert_name, 'w', bitOrder);
            end

            if metadata.numWrittenFreqs == numBands; continue; end   % if all freqs already written, move to next channel (not sure if this is really necessary)
            for iFreq = metadata.numWrittenFreqs + 1 : numBands
                [n,fo,mo,w] = firpmord(freqBands(iFreq, :), ...
                                       M, ...
                                       [1 1 1] * 0.01, ...
                                       Fs);
                b = firpm(n,fo,mo,w);
                metadata.b{iFreq} = b;
                
                disp(['Filtering ' header.channel(lfp_idx).name ...
                      ', wire ' num2str(repWire) ...
                      ', centered at ' num2str(center_freq(iFreq)) ' Hz']);
                
                tic
                lfpFilt = filtfilt(b, 1, volt_lfp);
                toc
                
                % now calculate the analytic signal
                tic
                lfp_analytic = hilbert(lfpFilt)';
                toc
                
                sigmax = max(abs([real(lfp_analytic); imag(lfp_analytic)]));
                metadata.sigmax(iFreq) = sigmax;
                
                scaled_analytic = (32767 / sigmax) * lfp_analytic;
                real_sig = int16(real(scaled_analytic));
                imag_sig = int16(imag(scaled_analytic));
                
%                 fwrite(fid, [real(lfp_analytic), imag(lfp_analytic)], 'double');
                fwrite(fid, [real_sig, imag_sig], 'int16');
                
                metadata.numWrittenFreqs = iFreq;
                save(metadata_filename, 'metadata');
                
            end
            sprintf('%f seconds to calculate channel', toc(chtic))
            
            fclose(fid);

            metadata.numWrittenChannels = iCh;
            save(metadata_filename, 'metadata');
            metadata.numWrittenFreqs = 0;
        end
        
    end
    
end
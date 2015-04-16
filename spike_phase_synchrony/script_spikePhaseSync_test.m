lfp_root  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal LFPs';
eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};
numEvents = length(eventList);
sta_width = 1;    % in seconds
winDuration = 0.1;  % in seconds, width
winStepSize = 0.05;  % advance window by 50 ms each time

[chDB_list, chDB_fnames] = get_chStructs_for_analysis;

implantID = implantID_from_ratID(chDB_list{1}(1:3));    % D20/IM166
channels = D20Ch_beta;   %eval( chDB_info.name );
spike    = D20spike;

lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
cd(lfp_directory);

iSpike = 60;
sp = spike{iSpike};

cp = initChanParams();
cp.channelName = spike.name(1:end-1);

chList = extractChannels(cp, channels);
ch     = channels(chList);

sessionDate = ch.date;

lfp_fileinfo = dir(['*_' sessionDate '.hsdf']);

%         if isempty(lfp_fileinfo); continue; end
        
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

all_fp = readHSD( lfp_fileName, ...
                  numRecordingSites, ...
                  header.dataOffset, ...
                  Fs, ...
                  [0, lfpDuration] );
           
repWire = getRepWise( ch );
lfp_idx = find( lfp_wireNums == repWire );
lfp.lfp = all_fp(lfp_idx, :)';
lfp.Fs  = Fs;

spike_ts = sp.timestamps.spike;

et_lfps = getEventTriggeredLFP(spike_ts, lfp, 
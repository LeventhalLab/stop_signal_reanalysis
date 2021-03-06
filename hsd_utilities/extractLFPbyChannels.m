function extractLFPbyChannels(channels, varargin)

targetFs = 500;
filtOrder = 1000;
dataType = 'int16';
saveDir = '';
openDirStart = '';
HSDblockSize = 100000;
overwriteOldLFP = 0;


% function to extract LFPs for a channel structure using the requested
% parameters (default target sampling rate 500 Hz)

for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'targetfs',
            targetFs = varargin{iarg + 1};
        case 'savedir',
            saveDir = varargin{iarg + 1};
        case 'opendir',
            openDirStart = varargin{iarg + 1};
        case 'filterorder',
            filtOrder = varargin{iarg + 1};
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'blocksize',
            HSDblockSize = varargin{iarg + 1};
        case 'overwritelfp',
            overwriteOldLFP = varargin{iarg + 1};
    end
end

bytes_per_sample = getBytesPerSample(dataType);
openDir = openDirStart;

% set channel parameters to extract any channel
% cp.task = -1;
% cp.locationName = 'any';
% cp.locationSubClass = 'any';
% cp.subject = 'any';
% cp.date = 'any';
% cp.tetrode = 'any';
% cp.session = 'any';
% cp.channelName = 'any';
cp = initChanParams();

% don't bother to calculate LFPs for reference ECoG or EMGs
cp.locationSubClass = {'emg', 'eeglam', 'ref'};
channels = excludeChannels(cp, channels);
cp.locationSubClass = 'any';

sessionList = getSessionsfromChannelDB(channels);
numSessions = length(sessionList);

for iSession = 1 : numSessions
    cp.session = sessionList{iSession}
    
    chList = extractChannels(cp, channels);
    c = channels(chList);
    numLFPchannels = length(chList);
    
    hsdName = c{1}.files.highSpeedData.file;
    if strcmp(hsdName(1),'\')
        hsdName = PCfn2macfn(hsdName);
    end
%     [pn, HSDfn, ext, ~] = fileparts(hsdName);
    [pn, HSDfn, ext] = fileparts(hsdName);
    if isempty(openDirStart)
        openDir = pn;
    else
        openDir = openDirStart;
    end
    
    sessionDir = fullfile(openDir, [HSDfn(1:14), 'a']);
    if ~exist(fullfile(sessionDir, [HSDfn ext]), 'file')
        disp(['File ' fullfile(sessionDir, HSDfn) ' does not exist.'])
        continue;
    end   
%     if ~exist(fullfile(openDir, [HSDfn ext]), 'file')
%         disp(['File ' fullfile(openDir, HSDfn) ' does not exist.'])
%         continue;
%     end   
    
    LFPfn = [HSDfn ext 'f'];
    
    if ~isempty(saveDir)
        pn = saveDir;
    elseif ~isempty(c{1}.files.lfp)
%         [pn, ~, ~, ~] = fileparts(c{1}.files.lfp.file);
        [pn, ~, ~] = fileparts(c{1}.files.lfp.file);
    end
    
    if exist(fullfile(pn, LFPfn), 'file') && ~overwriteOldLFP
        disp(['File ' fullfile(pn, LFPfn) ' already exists.']);
        continue;
    end
    
    hsdName = fullfile(openDir, [HSDfn ext]);
    hsdHeader = getHSDHeader(hsdName);
    hsdFs = hsdHeader.main.sampling_rate;
    dataOffset = hsdHeader.dataOffset;
    numHSDchannels = hsdHeader.main.num_channels;
    hsdInfo = dir(hsdName);
    numHSDsamples = (hsdInfo.bytes - dataOffset) / (numHSDchannels * bytes_per_sample);

    r = round(hsdFs / targetFs);
    lfp_Fs = hsdFs / r;
    overlapSize = ceil(filtOrder / r) * r;
    
    lfpHeader = hsdHeader;
    lfpHeader.main.downsampled_rate = double(lfp_Fs);
%     lfpHeader.comment = ['HSD downsampled by ' num2str(r) ' from ' num2str(hsdFs)];
    lfpHeader.dataFilename = LFPfn;
    
    % the line below makes sure that there are only as many channels
    % contained in the lfpfHeader.channels structure as there are LFPs that
    % will be made
    lfpHeader.channel = hsdHeader.channel(1);
    chIdx = 0;
    for iCh = 1 : numLFPchannels
        tempRepWire = getRepWire(c{iCh});
        if tempRepWire == 0
            % no good wires for this channel
            continue;
        end
        chIdx = chIdx + 1;
        repWire(chIdx) = tempRepWire;
        lfpHeader.channel(chIdx) = hsdHeader.channel(tempRepWire);
    end
    numLFPchannels = chIdx;
    lfpHeader.main.num_channels = uint16(numLFPchannels);
    
    LFPfn = fullfile(pn, LFPfn);
    
    lfpHeader.main.version = 3;
    writeHSDheader(LFPfn, lfpHeader);
   
    numLFPsamples = ceil(numHSDsamples / r);
    lfpData = zeros(numLFPchannels, numLFPsamples);
    
    HSDblockSize = ceil(HSDblockSize / r) * r;
    numBlocks = ceil(numHSDsamples / HSDblockSize);
    LFPblockSize = HSDblockSize / r;
    HSD = zeros(numHSDchannels, HSDblockSize);
    
    % do the first block separate from the rest, since there will be some
    % overlap for the rest of the blocks
    fin = fopen(hsdName, 'r');
    fseek(fin, dataOffset, 'bof');
    
    LFPstart = 1;
    LFPend = LFPstart + LFPblockSize - 1;
    % do the first block
    disp(['Block 1 of ' num2str(numBlocks)]);
    HSD = fread(fin, [numHSDchannels, HSDblockSize], 'int16', 0, 'b');
%     activeHSD = HSD(repWire,:);
    currentLFP = zeros(numLFPchannels, LFPblockSize);
    for iCh = 1 : numLFPchannels
        currentLFP(iCh, :) = ...
            decimate(HSD(repWire(iCh), :), r, filtOrder, 'fir');
    end
    lfpData(:, LFPstart:LFPend) = currentLFP;
    
    LFPstart = LFPend + 1;
    HSDblockSize = HSDblockSize + overlapSize;
    activeHSD = zeros(numLFPchannels, HSDblockSize);
    currentLFP = zeros(numLFPchannels, HSDblockSize / r);
    LFPblockStart = overlapSize / r + 1;
    tic
    for iBlock = 2 : numBlocks
        disp(['Block ' num2str(iBlock) ' of ' num2str(numBlocks) '; session # ' num2str(iSession) ' of ' num2str(numSessions) ' (' sessionList{iSession} ')']);
        rewindSize = -overlapSize * numHSDchannels * bytes_per_sample;
        fseek(fin, rewindSize, 'cof');
        
        HSD = fread(fin, [numHSDchannels, HSDblockSize], 'int16', 0, 'b');
%         activeHSD = HSD(repWire, :);
        if iBlock < numBlocks
            LFPend = LFPstart + LFPblockSize - 1;
        else
            LFPend = size(lfpData, 2);
        end
        
        if iBlock == numBlocks
            clear currentLFP;
        end
        for iCh = 1 : numLFPchannels
            currentLFP(iCh,:) = ...
                decimate(HSD(repWire(iCh), :), r, filtOrder, 'fir');
        end
        
        lfpData(:, LFPstart:LFPend) = currentLFP(:, LFPblockStart : size(currentLFP, 2));
        LFPstart = LFPend + 1;
        
    end
    fclose(fin);
    
    lfpData = int16(lfpData);
    
    fout = fopen( LFPfn, 'r+');
    fseek(fout, dataOffset, 'bof');
    fwrite(fout, lfpData, 'int16', 0, 'b');
    fclose(fout);
    toc
    
end
    
    
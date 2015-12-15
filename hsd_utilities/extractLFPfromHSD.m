function extractLFPfromHSD(hsdName, varargin)

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

if strcmp(hsdName(1),'\')
    hsdName = PCfn2macfn(hsdName);
end
[pn, HSDfn, ext, ~] = fileparts(hsdName);
if isempty(openDirStart)
    openDir = pn;
end
if ~exist(fullfile(openDir, [HSDfn ext]), 'file')
    disp(['File ' fullfile(openDir, HSDfn) ' does not exist.'])
    return;
end    

LFPfn = [HSDfn ext 'f'];

if ~isempty(saveDir)
    pn = saveDir;
end

if exist(fullfile(pn, LFPfn), 'file') && ~overwriteOldLFP
    disp(['File ' fullfile(pn, LFPfn) ' already exists.']);
    return;
end

hsdHeader = getHSDHeader(hsdName);
hsdFs = hsdHeader.main.sampling_rate;
dataOffset = hsdHeader.dataOffset;
numHSDchannels = hsdHeader.main.num_channels;
hsdInfo = dir(hsdName);
numHSDsamples = (hsdInfo.bytes - dataOffset) / (numHSDchannels * bytes_per_sample);

r = round(hsdFs / targetFs);
lfpFs = hsdFs / r;
overlapSize = ceil(filtOrder / r) * r;

lfpHeader = hsdHeader;
lfpHeader.main.sampling_rate = uint16(round(lfpFs));
lfpHeader.comment = ['HSD downsampled by ' num2str(r) ' from ' num2str(hsdFs)];
lfpHeader.dataFilename = LFPfn;

numLFPchannels = length(lfpHeader.channel);

LFPfn = fullfile(pn, LFPfn);
writeHSDheader(LFPfn, lfpHeader);

numLFPsamples = ceil(numHSDsamples / r);
lfpData = zeros(numLFPchannels, numLFPsamples);

HSDblockSize = ceil(HSDblockSize / r) * r;
numBlocks = ceil(numHSDsamples / HSDblockSize);
LFPblockSize = HSDblockSize / r;

% do the first block separate from the rest, since there will be some
% overlap for the rest of the blocks
fin = fopen(hsdName, 'r');
fseek(fin, dataOffset, 'bof');

LFPstart = 1;
LFPend = LFPstart + LFPblockSize - 1;
% do the first block
disp(['Block 1 of ' num2str(numBlocks)]);
HSD = fread(fin, [numHSDchannels, HSDblockSize], 'int16', 0, 'b');

currentLFP = zeros(numLFPchannels, LFPblockSize);
parfor iCh = 1 : numLFPchannels
    currentLFP(iCh, :) = ...
        decimate(HSD(iCh, :), r, filtOrder, 'fir');
end
lfpData(:, LFPstart:LFPend) = currentLFP;

LFPstart = LFPend + 1;
HSDblockSize = HSDblockSize + overlapSize;

currentLFP = zeros(numLFPchannels, HSDblockSize / r);
LFPblockStart = overlapSize / r + 1;
tic
for iBlock = 2 : numBlocks
    disp(['Block ' num2str(iBlock) ' of ' num2str(numBlocks)]);
    rewindSize = -overlapSize * numHSDchannels * bytes_per_sample;
    fseek(fin, rewindSize, 'cof');

    HSD = fread(fin, [numHSDchannels, HSDblockSize], 'int16', 0, 'b');

    if iBlock < numBlocks
        LFPend = LFPstart + LFPblockSize - 1;
    else
        LFPend = size(lfpData, 2);
    end

    if iBlock == numBlocks
        clear currentLFP;
    end
    parfor iCh = 1 : numLFPchannels
        currentLFP(iCh,:) = ...
            decimate(HSD(iCh, :), r, filtOrder, 'fir');
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
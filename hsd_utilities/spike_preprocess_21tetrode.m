function spike_preprocess_21tetrode(varargin)
%
% usage: spike_preprocess_21tetrode(varargin)
%
% This function takes a .hsd file from a 21-tetrode implant and will
% wavelet-filter the data in preparation for spike sorting. It will write a
% separate file for each tetrode/ref (stereotrode) labeled with the
% original .hsd name, with an added "_Txx" or "_Rxx" indicating the tetrode
% or ref #. Each file will contain an appropriate number of wires (4 for
% tetrodes, 2 for refs) written big-endian, and should be readable with
% standard readHSD commands.
%
% INPUTS:
%
% VARARGINS:
%       'hsdname' - name of the .hsd file from which to get raw data
%       'targetdir' - name of directory in which to store .hsdw files
%       'machineformat' - 'b' (big-endian) or 'l' (little-endian). Labview
%           on Windows records big-endian; Sunjay's human data are recorded
%           little-endian, at least as of 11/18/2011
%       'datatype' - default int16
%       'fs' - sampling rate; over-ridden by header if available
%       'tetrode' - array of tetrode numbers on which to perform wavelet
%           filtering. Any tetrodes not listed in header are skipped 
%       'ref' - array of ref (stereotrode) numbers to wavelet filter.
%
% OUTPUTS:

blockSize        = 100000;
bufferLength     = 10000;   % to prevent edge effects
netBlockSize     = blockSize + 2 * bufferLength;
dataType         = 'int16';
bytes_per_sample = getBytesPerSample(dataType);

fn = '';
targetDir = '';
machineFormat = 'b';
dataType = 'int16';
tetrodes = [];
refs = [];   % stereotrodes

for iarg = 1 : 2 : nargin
    
    switch lower(varargin{iarg})
        case 'hsdname',
            fn = varargin{iarg + 1};
        case 'targetdir',
            targetDir = varargin{iarg + 1};
        case 'machineformat',
            machineFormat = varargin{iarg + 1};
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'tetrode',
            tetrodes = varargin{iarg + 1};
        case 'ref',
            refs = varargin{iarg + 1};

    end
    
end

if isempty(fn)
    [fname, pn, ~] = uigetfile('*.hsd');
    fn = fullfile(pn, fname);
else
    [pn, fname, ~, ~] = fileparts(fn);
end
if isempty(targetDir)
    targetDir = uigetdir(pn);
end

if ~targetDir
    disp('no target directory');
    return;
end

if ~isdir(targetDir)
    mkdir(targetDir);
end

hasHeader = checkHSDheader( fn );
if ~hasHeader
    disp(['no header for ' fn]);
    return;
end

header = getHSDHeader( fn );
ch = header.channel;
numCh = length(ch);
chTypes = [ch.channel_type];           % list of channel types for each wire
                                                         %   1- tetrode
                                                         %   2- ref (stereotrode)
                                                         %   3- eeg/emg/true reference
chNums = [ch.channel_number];      % list of channel numbers (ie, tetrode 7, EEG 2, etc.)
dataOffset = header.dataOffset;

% if specific tetrodes haven't been selected, use all the tetrodes
if isempty(tetrodes)
    tetrodes = unique(chNums(chTypes == 1));
end
% if specific refs haven't been selected, use all the refs (stereotrodes)
if isempty(refs)
    refs = unique(chNums(chTypes == 2));
end
if strcmpi(refs, 'none')
    refs = [];
end
if strcmpi(tetrodes,'none')
    tetrodes = [];
end

validTypes = (chTypes == 1);  % find all tetrodes
numTetrodes = length(tetrodes);
numRefs = length(refs);
fout = zeros(1, numTetrodes + numRefs);
wires = cell(1, numTetrodes + numRefs);
hsdwHeaders = cell(1, numTetrodes + numRefs);
saveName = cell(1, numTetrodes + numRefs);
for iTet = 1 : numTetrodes
    wires{iTet} = find(validTypes & (chNums == tetrodes(iTet)));
    if ~isempty(wires{iTet})
        tetName = ['T' sprintf('%02d', tetrodes(iTet))];
        saveName{iTet} = fullfile(targetDir, [fname '_' tetName '.hsdw']);
        hsdwHeaders{iTet} = createHSDWHeader(header, 1, tetrodes(iTet));
        writeHSDheader(saveName{iTet}, hsdwHeaders{iTet});
        fout(iTet) = fopen(saveName{iTet}, 'r+');
        fseek(fout(iTet), dataOffset, 'bof');
    end
end

validTypes = (chTypes == 2);  % find all stereotrodes
for iRef = 1 : numRefs
    wires{numTetrodes + iRef} = find(validTypes & (chNums == refs(iRef)));
    if ~isempty(wires{numTetrodes + iRef})
        refName = ['R' sprintf('%02d', refs(iRef))];
        saveName{numTetrodes + iRef} = fullfile(targetDir, [fname '_' refName '.hsdw']);
        hsdwHeaders{numTetrodes + iRef} = createHSDWHeader(header, 2, refs(iRef));
        writeHSDheader(saveName{numTetrodes + iRef}, hsdwHeaders{numTetrodes + iRef});
        fout(numTetrodes + iRef) = fopen(saveName{numTetrodes + iRef}, 'r+');
        fseek(fout(numTetrodes + iRef), dataOffset, 'bof');
    end
end

totChannels = length(wires);
% now load in the .hsd data in blocks and wavelet filter it
fileInfo   = dir(fn);
numBytes   = fileInfo.bytes;

data_bytes = numBytes - dataOffset;
hsdSamples = data_bytes / (bytes_per_sample * numCh);

numBlocks  = ceil(hsdSamples / blockSize);

fin = fopen(fn, 'r');
fseek(fin,dataOffset,'bof');
for iBlock = 1 : numBlocks
    disp(['wavelet filtering ' fname ', block ' num2str(iBlock) ' of ' num2str(numBlocks)]);
    
    
    hsd = fread(fin, [numCh, netBlockSize], dataType, 0, machineFormat);
    
    for iChannel = 1 : totChannels
        if isempty(wires{iChannel}); continue; end
        
        hsd_to_filt = hsd(wires{iChannel}, :);
        fdata = wavefilter(hsd_to_filt, 6);
        
        if iBlock == 1
            fdata = fdata(:, 1:blockSize);
        elseif iBlock == numBlocks
            fdata = fdata(:, bufferLength+1:end);
        else
            fdata = fdata(:, bufferLength+1:bufferLength+blockSize);
        end
        
        fwrite(fout(iChannel), fdata, dataType, 0, 'b');
    
    end    % for iChannel...

    if iBlock == 1
        % if first block, rewind to bufferLength samples before the end of
        % the previous block. This is because 2*bufferLength extra samples
        % are read on the first pass, but bufferLength samples are read in
        % at either side of subsequent passes
        rewindLength = -(bufferLength * 3) * numCh * bytes_per_sample;
    else
        % if not the first block, rewind to bufferLength samples before the
        % end of the previous block
        rewindLength = -(bufferLength * 2)* numCh * bytes_per_sample;
    end
    if iBlock < numBlocks
        fseek(fin, rewindLength, 'cof');
    end


    
end
fclose(fin);
for iFile = 1 : totChannels
    try
        fclose(fout(iFile));
    catch
    end
end

end    % function
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newHeader = createHSDWHeader( oldHeader, channelType, channelNum )
%
% INPUTS:
%   oldHeader - header from the original .hsd file
%   channelType - tetrode (1) or ref (2)
%   channelNum - channel number (ie, tetrode 7, ref 2, etc.). Note this is
%       confusing as sometimes a single wire is called a channel, and sometimes
%       the whole tetrode/ref is called a channel. What can you do? This is the
%       terminology I was stuck with. -DL
%
% OUTPUTS:
%   newHeader - new header; same as old header but with channel information
%       replaced with information specific to the tetrode/ref of interest
%

newHeader = oldHeader;

ch = oldHeader.channel;
chTypes = [ch.channel_type];           % list of channel types for each wire
                                                         %   1- tetrode
                                                         %   2- ref (stereotrode)
                                                         %   3- eeg/emg/true reference
chNums = [ch.channel_number];      % list of channel numbers (ie, tetrode 7, EEG 2, etc.)

validTypes = (chTypes == channelType);  % find all tetrodes or all refs
wireList = (validTypes & (chNums == channelNum));

newHeader.main.num_channels = length(find(wireList));
newHeader.channel = ch(wireList);

end    % function newHeader...

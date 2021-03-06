function [out] = readSingleWire( fn, wireNum, varargin )
%
% usage: [out] = readSingleWire( fn, wireNum, varargin )
%
% function to read in a single wire of raw data from an hsd/hsdw file given
% the filename, number of channels, data offset, sampling rate, and
% timelimits
%
% USAGE:
% fn = 'myWonderfulData.hsd';
% numChannels = 81;   (for a 21-tetrode drive)
% dataOffset = 20*1024; % this accounts for the header data
%
% varargs:
%   'timelimits' - 2-element vector with start time and end time to read in
%      (in seconds). Default is the full recording.
%   'datatype' - type of data; default int16
%   'fs' - sampling rate in Hz; if not provided, the sampling rate is read
%      from the .hsd header, or calculated from the .hsdf (LFP) header

dataType = 'int16';
timeLimits = [0 0];
Fs = 0;
numSampsGiven = 0;
bitOrder = 'b';

for iarg = 1 : 2 : nargin - 2
    
    switch lower(varargin{iarg})
        
        case 'timelimits',
            timeLimits = varargin{iarg + 1};
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'fs',
            Fs = varargin{iarg + 1};
        case 'starttime_numsamps',
            timeLimits(1) = varargin{iarg + 1}(1);
            timeLimits(2) = 1;
            numSampsGiven = 1;
            numSamples    = varargin{iarg + 1}(2);
        case 'bitorder',
            bitOrder = varargin{iarg + 1};
    end
    
end

if ~timeLimits(2)   % default time limits are the entire file
    timeLimits = [0 getHSDlength( 'filename', fn, 'datatype', dataType )];
end

bytes_per_sample = getBytesPerSample( dataType );
hsd_header = getHSDHeader( fn );
if ~Fs
    % sampling rate not specified in input to this function
    [~, ~, ext, ~] = fileparts(fn);
    if strcmpi(ext,'.hsdf')
        % this is a .hsdf (LFP) file
        if isfield(hsd_header.main, 'downsampled_rate')
            if hsd_header.main.downsampled_rate < 100
                Fs = lfpFs(hsd_header);
            else
                Fs = hsd_header.main.downsampled_rate;
            end
        else
            Fs = lfpFs(hsd_header);
        end
    else
        % this is a raw data (.hsd) file
        Fs = hsd_header.main.sampling_rate;
    end
end
numChannels = hsd_header.main.num_channels;
dataOffset = hsd_header.dataOffset;

% check which wire is being pulled out (for example, the nth row in the
% data array may not be the nth wire)
channelNum = 0;
for iChannel = 1 : numChannels
    if hsd_header.channel(iChannel).original_number == wireNum
        channelNum = iChannel;
        break;
    end
end

if ~ channelNum
    % there were no matches for the target wire number in the hsd header
    disp(['Error in readSingleWire - wire ' num2str(wireNum) ' was not found in the list of channels in the file header.']);
    out = 0;
    return;
end

fin = fopen(fn, 'r');

startSample = floor( Fs * timeLimits(1));
if ~numSampsGiven
    numSamples = ceil(Fs * range(timeLimits));
end

startPosition = dataOffset + (startSample * numChannels + (channelNum - 1)) * bytes_per_sample;
skipBytes = (numChannels - 1) * bytes_per_sample;

fseek( fin, startPosition, 'bof' );
% if strcmp(fn(length(fn) - 3 : length(fn)), '.hsd')
    [out, number_elements_read] = fread(fin, numSamples, dataType, skipBytes, bitOrder);
% elseif strcmp(fn(length(fn) - 4 : length(fn)), '.hsdf')
%     [out, number_elements_read] = fread(fin, numSamples, dataType, skipBytes);
% end

fclose(fin);
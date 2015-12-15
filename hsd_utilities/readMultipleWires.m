function [data] = readMultipleWires(fn, wireNums, varargin)

% returns data in an m x n array, where each column is one wire's worth of
% data

dataType = 'int16';
timeLimits = [0 0];
Fs = 0;

for iarg = 1 : 2 : nargin - 2
    
    switch lower(varargin{iarg})
        
        case 'timelimits',
            timeLimits = varargin{iarg + 1};
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'fs',
            Fs = varargin{iarg + 1};
    end
    
end

numWires = length(wireNums);

if ~timeLimits(2)   % default time limits are the entire file
    timeLimits = [0 getHSDlength( 'filename', fn, 'datatype', dataType )];
end

hsd_header = getHSDHeader( fn );
if ~Fs
    % sampling rate not specified in input to this function
    [~, ~, ext, ~] = fileparts(fn);
    if strcmpi(ext,'.hsdf')
        % this is a .hsdf (LFP) file
        Fs = lfpFs(hsd_header);
    else
        % this is a raw data (.hsd) file
        Fs = hsd_header.main.sampling_rate;
    end
end

numSamples = ceil(Fs * range(timeLimits));
data = zeros(numSamples, numWires);

for iCh = 1 : numWires
    if wireNums(iCh) == 0
        data(:, iCh) = zeros(numSamples, 1);
    else
        data(:, iCh) = readSingleWire(fn, wireNums(iCh), ...
                   'timelimits', timeLimits, 'datatype', dataType);
    end
end
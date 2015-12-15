function rawData = readHSDSamples( fn, numWires, dataOffset, startSample, numSamples, varargin )
% function to read in raw data from an hsd/hsdw file given the filename,
% number of channels, data offset, sampling rate, starting sample, and
% number of samples to read. This differs from readHSD, which reads in data
% based on a time range. This function is used to prevent mismatches in
% array sizes due to rounding errors in the number of samples to read
%
% data = readHSD( filename, numWires, dataOffset, Fs, timeLimits )
%
% USAGE:
% filename = 'myWonderfulData.hsd';
% numWires = 81;
% dataOffset = 20*1024; % this accounts for the header data
% Fs = 31250;
% timeLimits = [0 30]; % read the first 30 seconds
%

bytes_per_sample = 2;

for iarg = 1 : 2 : nargin - 5
    
    switch lower(varargin{iarg})
        
        case 'bytespersample',
            bytes_per_sample = varargin{iarg + 1};
            
        case 'datatype',
            dataType = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...

switch bytes_per_sample
    case 2,
        dataType = 'int16';
end    % end switch bytes_per_sample

fin = fopen(fn, 'r');

startPosition = dataOffset + startSample * numWires * bytes_per_sample;
fseek(fin, startPosition, 'bof');

if fn(length(fn)) == 'd'    % this is a .hsd file
    
    if or(strcmp(fn(length(fn) - 7 : length(fn)), 'wire.hsd'), ...
          strcmp(fn(length(fn) - 7 : length(fn)), '.lfp.hsd'));    % this is a single wire file or lfp file
        [rawData, number_elements_read] = fread(fin, [numWires, numSamples], dataType);
    else            % this is a .hsd or .hsdf file
        [rawData, number_elements_read] = fread(fin, [numWires, numSamples], dataType, 0,'b');
    end
    
else    % this is a .hsdw file
    
    [rawData, number_elements_read] = fread(fin, [numWires, numSamples], dataType,0,'b');
    
end    
    

fclose(fin);
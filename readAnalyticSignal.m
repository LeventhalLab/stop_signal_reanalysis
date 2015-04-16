function analyticSignals = readAnalyticSignal( fn, metadata, timelimits, freqIdx )
%
% function to read analytic signals written into binary files as double
% precision complex numbers
%
% usage: 
%
% INPUTS:
%   fn - filename (including path) of the binary file containing the
%        analytic signals
%   metadata - metadata structure stored with the analytic signals
%   timelimits - 2 element vector containing the start and end time for
%       which to pull down the analytic signal
%   freqIdx - indices of the frequency bands to pull down; to figure out
%       which frequency bands correspond with which index, look at the
%       metadata.freqBands array, which should be an m x 4 array where m is
%       the number of different frequency bands and the rows contain the
%       frequency ranges used for making the filters
%
% OUTPUTS:
%   analyticSignals - 

bytes_per_sample = 2;    % 16-bit integers

if ~exist(fn, 'file')
    error([fn ' could not be found.']);
end

Fs          = metadata.Fs;
startSamp   = round(timelimits(1) * Fs);
sampsToRead = round(range(timelimits) * Fs);
totalSamps  = metadata.numSamps;
numFreqs    = length(freqIdx);

% figure out the index of the current channel
% fn_analytic_idx   = strfind(fn, 'analytic_');
% fn_chname_end_idx = strfind(fn, '.bin');
% chname_from_fn  = fn(fn_analytic_idx + 9 : fn_chname_end_idx-1);
% chIdx = find(strcmp(chname_from_fn, metadata.chNames));

re_startBytes = ((freqIdx-1) * 2 * totalSamps + startSamp) * bytes_per_sample;
im_startBytes = (((freqIdx-1) * 2 + 1) * totalSamps + startSamp) * bytes_per_sample;

re_analytic = double(zeros(sampsToRead, numFreqs));
im_analytic  = double(zeros(sampsToRead, numFreqs));

fid = fopen(fn, 'r');

for iFreq = 1 : numFreqs

    fseek(fid, re_startBytes(iFreq), 'bof');
    re_analytic(:, iFreq) = double(fread(fid, sampsToRead, 'int16', 0, metadata.bitOrder));
    re_analytic(:, iFreq) = re_analytic(:, iFreq) * metadata.sigmax(freqIdx(iFreq)) / 32767;
    
    fseek(fid, im_startBytes(iFreq), 'bof');
    im_analytic(:, iFreq) = double(fread(fid, sampsToRead, 'int16', 0, metadata.bitOrder));
    im_analytic(:, iFreq) = im_analytic(:, iFreq) * metadata.sigmax(freqIdx(iFreq)) / 32767;
end

fclose(fid);

analyticSignals = complex(re_analytic, im_analytic);
function [MI, low_freqs, high_freqs] = phase_amp_comodugram( ch, varargin )
%
% 
%
% usage:
%
% INPUTS:
%   ch - 
%
% OUPUTS:
%

low_freq_range  = [0, 20];
high_freq_range = [8, 100];

phasebins = 10:20:360;

bitOrder = 'b';

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
lfp_root          = '/Volumes/Recordings/dan/Leventhal Neuron 2012_summary/Leventhal Neuron 2012_LFPs';
hilbert_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';

for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'lowfreqrange',
            low_freq_range = varargin{iarg + 1};
        case 'highfreqrange',
            high_freq_range = varargin{iarg + 1};
        case 'phasebins',
            phasebins = varargin{iarg + 1};
    end
end
            
% navigate to the directory in which the analytic signals are stored
implantID = implantID_from_ratID(ch.name(1:3));
subject_hilbertDir = fullfile(hilbert_directory, [implantID '_hilbert']);

hilbert_sessionDir = fullfile(subject_hilbertDir, ch.session);

if ~exist(hilbert_sessionDir, 'dir')
    error([hilbert_sessionDir ' could not be found.']);
end

metadata_filename = [ch.session 'hilbert_metadata.mat'];
metadata_filename = fullfile(hilbert_sessionDir, metadata_filename);

if ~exist(metadata_filename, 'file')
    error('metadata file not found.');
end

load(metadata_filename);

hilbert_name = ['analytic_' ch.name '.bin'];
hilbert_name = fullfile(hilbert_sessionDir, hilbert_name);

centerFreqs = mean(metadata.freqBands, 2);
low_freq_idx  = find(centerFreqs > low_freq_range(1) & centerFreqs < low_freq_range(2));
high_freq_idx = find(centerFreqs > high_freq_range(1) & centerFreqs < high_freq_range(2));

low_freqs  = centerFreqs(low_freq_idx);
high_freqs = centerFreqs(high_freq_idx);

num_low_freq  = length(low_freq_idx);
num_high_freq = length(high_freq_idx);

MI = zeros(num_low_freq, num_high_freq);
U = ones(1, length(phasebins)) / length(phasebins);

for if1 = 1 : num_low_freq
    disp([ch.name ', low freqency = ' num2str(mean(metadata.freqBands(if1,:),2))]);
    
    as1 = readAnalyticSignal(hilbert_name, metadata, ...
                            [0 metadata.duration], ...
                            low_freq_idx(if1));
    for if2 = 1 : num_high_freq
        
        % no point in looking for nesting of a lower frequency oscillations
        % inside a higher frequency one
        if low_freq_idx(if1) >= high_freq_idx(if2); continue; end
        
%         disp('loading analytic signals...')
%         tic
        as2 = readAnalyticSignal(hilbert_name, metadata, ...
                                 [0 metadata.duration], ...
                                 high_freq_idx(if2));
%         toc
        
%         disp('calculating phase-amplitude histogram...')
%         tic
        [y, phasebins] = phase_amp_histogram( as1, as2 );
%         toc
        y = y ./ sum(y);    % normalize over sum of mean amplitudes across bins
        MI(if1, if2) = KLdistance(y, U) / log10(length(phasebins));
        
    end
end
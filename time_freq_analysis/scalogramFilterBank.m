function [filterBank] = scalogramFilterBank(f, Fs, numSamples, varargin)
%
% this version is slightly different from Alex's original in that phase
% information is retained (the output is complex)
%
% Usage:
% W = calculateGabor_EnMasse(data, 'sigma', 0.2, 'fpass', [1 100], 'numfreqs',
%                   100, 'Fs', 1024)
%
% Input:
% f - vector containing center frequencies at which to calculate Morlet
%     wavelets
% Fs - sampling rate
% numSamples - number of samples in the original data stream
%
% Output:
% filterBank - a 2D array where each column is the fourier transform of a
%              complex Gabor wavelet with a center frequency given by the
%              frequencies in f
%
% doplot: do a plot of the 
% showwavelet: show the gabor wavelet (kernel, function, whatever) used for
%         analysis. plots the wavelet using a center frequency in the 
%         middle of the fpass range. (off by default, I use this to get a 
%         sense if my wavelet looks okay)
%
% Description: 
% Create a filter bank for decomposing a time signal (or set of time
% signals) into time-frequency components using Gabor (Morlet) wavelets

doplot = false;

samples_to_pad = 0;   % number of samples to pad the filter bank with on each side of the center Gabor wavelet
if size(varargin)>0
    for iarg= 1:2:length(varargin),   % assume an even number of varargs
        switch lower(varargin{iarg}),
			case {'doplot', 'plot'}
				doplot = varargin{iarg+1};
            case {'numpadsamples'}
                samples_to_pad = varargin{iarg+1};
        end % end of switch
    end % end of for iarg
end

% when we actually calculate the transform, we will pad the data array with
% zeros to prevent edge effects.
if samples_to_pad == 0
    numPaddedSamples = numSamples + 2 * round(numSamples/2);
else
    numPaddedSamples = numSamples + 2 * samples_to_pad;
end

%% Create the filter bank (where each kernel is the length of the signal input)
t = linspace(-numPaddedSamples/2, numPaddedSamples/2, numPaddedSamples)/Fs;

recip_gsigma = f ./ 0.849;
gaussWindow = exp( -0.5*(t'*recip_gsigma).^2 );
sinusoid_matrix = exp(1i*2*pi*t'*f);
t_filterBank = (1/pi^0.25) * gaussWindow .* sinusoid_matrix;

filterBank = fft(t_filterBank); % move the kernels into frequency space

%%
if doplot
   % add code here to plot the filterbank if requested
end

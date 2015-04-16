function [y, varargout] = phase_amp_histogram( fp, fa, varargin )
%
% function to compute the phase-amplitude coupling histogram for two
% complex time-series f1 and f2
%
% usage:
%
% INPUTS:
%   fp - time-series from which to extract phase
%   fa - time-series from which to extract amplitude
%
% VARARGINS:
%   bins - centers of phase bins. Default 20 degree bins from 0 to 360
%   nbins - number of evenly spaced bins to create
%   angleunit - angular units to use. Default "degrees", alternative
%       "radians"
%
% OUTPUTS:
%   y - the phase-amplitude histogram
%
% VARARGOUTS:
%   bins - phase bin centers

angleunit = 'degrees';
phasebins = 10 : 20 : 350;
nbins     = length(phasebins);

user_supplied_bins = false;

for iarg = 1 : 2 : nargin - 2
    switch lower(varargin{iarg})
        case 'angleunit',
            angleunit = varargin{iarg + 1};
        case 'bins',
            phasebins = varargin{iarg + 1};
            user_supplied_bins = true;
        case 'nbins',
            phasebins = 360/varargin{iarg+1}/2 : ...
                        360/varargin{iarg+1} : ...
                        360;
    end
    
end

nbins = length(phasebins);
% convert bins from degrees to radians, unless supplied by the user in
% radians
if ~user_supplied_bins || (user_supplied_bins && strcmpi(angleunit, 'degrees'))
    phasebins = phasebins * pi / 180;
end

% get the phase series
phase_angles = angle( fp );
phase_angles(phase_angles < 0) = phase_angles(phase_angles < 0) + 2 * pi;

% get the amplitude series
sig_amp = abs(fa);

bin_lim = zeros(1,2);
y = zeros(1, nbins);
for iBin = 1 : nbins
    if iBin == 1
        bin_lim(1) = 0;
    else
        bin_lim(1) = mean(phasebins( iBin-1 : iBin ));
    end
    if iBin == nbins
        bin_lim(2) = 2 * pi;
    else
        bin_lim(2) = mean(phasebins( iBin : iBin + 1));
    end
    
    y(iBin) = mean(sig_amp(phase_angles >= bin_lim(1) & phase_angles < bin_lim(2)));
    
end

varargout{1} = phasebins;

end
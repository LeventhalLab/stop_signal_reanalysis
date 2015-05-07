function [newFs, varargout] = lfpFs( lfp_header )
%
% function to pull out the true sampling rate from a .hsdf file header. For
% historical reasons, the downsampling ratio is stored in the comment in
% the phrase "HSD downsampled by xx from yyyyy", where xx is the
% downsampling ration (I've been using 30), and yyyyy is the original
% sampling rate (typically 31250).
%
% INPUTS:
%   lfp_header - the header from a .hsdf with the standard comment as above
%
% OUTPUTS:
%   newFs - the sampling rate for the LFP (the downsampled rate) as a
%       double
% varargouts:
%   r - the downsampling ratio
%   oldFs = the original sampling rate

if isfield(lfp_header.main, 'downsampled_rate')
    newFs = lfp_header.main.downsampled_rate;
end

if newFs > 100
    r = 0;
    oldFs = 0;
    return;
end
        
searchString = 'downsampled by ';
rLoc = findstr(lfp_header.comment, searchString) + length(searchString);
searchString = 'from ';
fLoc = findstr(lfp_header.comment, searchString) + length(searchString);

r = str2num(lfp_header.comment(rLoc : rLoc + 1));

oldFs = str2num(lfp_header.comment(fLoc : length(lfp_header.comment)));

newFs = oldFs / r;
varargout(1) = {r};
varargout(2) = {oldFs};
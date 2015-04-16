function D = KLdistance(P, Q)
%
% usage: D = KLdistance(P, Q)
%
% INPUTS:
%   P - probability distribution 1
%   Q - probability distribution 2 (must be same length as P)
%
% OUTPUTS:
%   D - Kullback-Leibler distance of P from Q

% make sure both are column vectors
if length(P) > size(P,1); P = P'; end
if length(Q) > size(Q,1); Q = Q'; end

D = nansum(P.*log10(P./Q));
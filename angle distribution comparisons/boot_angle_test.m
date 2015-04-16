function p = boot_angle_test(theta1, theta2, varargin)
%
% usage: 
%
% INPUTS:
%
% OUTPUTS:
%
% see Stark and Abeles, "Applying resampling methods to neurophysiological
%   data", J Neurosci Methods, 2005

nBoot = 1000;

for iarg = 1 : 2 : nargin - 2
    switch lower(varargin{iarg})
        case 'nboot',
            nBoot = varargin{iarg + 1};
    end
end

n1 = length(theta1); n2 = length(theta2);
if n1 > size(theta1, 1)
    theta1 = theta1';
end
if n2 > size(theta2, 1)
    theta2 = theta2';
end

cplx_theta1 = exp(1i*theta1);
cplx_theta2 = exp(1i*theta2);

fullDist = [theta1; theta2];
testStat = abs(mean(cplx_theta1) - mean(cplx_theta2));

bootStat = zeros(1, nBoot);
for iBoot = 1 : nBoot
    bootPerm = randperm(length(fullDist));
    
    cplx_resamp1 = exp(1i * fullDist(bootPerm(1:n1)));
    cplx_resamp2 = exp(1i * fullDist(bootPerm(n1+1:end)));
    bootStat(iBoot) = abs(mean(cplx_resamp1) - mean(cplx_resamp2));
    
end

p = 1 - (find(bootStat > testStat, 1, 'first') / nBoot);
if isempty(p); p = 1; end
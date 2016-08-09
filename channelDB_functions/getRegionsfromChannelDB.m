function [regionList] = getRegionsfromChannelDB( channels )
%
% function to extract the set of individual regions from a channel DB
% structure
%
% INPUTS:
%   channels - a channel structure
%
% VARARGINs:
%   onlyValidRegions - whether or not to include only regions recorded from
%      channels marked "valid"
%
% OUTPUTS:
%   regionList - a cell array containing the list of sessions
%      included in the "channels" channel DB


if isempty(channels)
    regionList = {};
    return;
end

regionList{1} = channels{1}.location.name;
numRegions = 1;
for iCh = 1 : length(channels)
    
    if max(strcmpi( channels{iCh}.location.name, regionList ))
        % this region has already been counted
        continue;
    end
    
    numRegions = numRegions + 1;
    regionList{numRegions} = channels{iCh}.location.name;
    
end
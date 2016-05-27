function [subclassList] = getSubclassesfromChannelDB( channels )
%
% usage: [sessionList] = getSubclassesfromChannelDB( channels )
%
% function to extract the set of individual subclasses from a channel DB
% structure
%
% INPUTS:
%   channels - a channel structure
%
% OUTPUTS:
%   subclassList - a cell array containing the list of location subclasses
%      included in the "channels" channel DB

subclassList{1} = channels{1}.location.subclass;
num_subClasses = 1;
for iCh = 1 : length(channels)
    
    if max(strcmpi( channels{iCh}.location.subclass, subclassList ))
        % this session has already been counted
        continue;
    end
    
    num_subClasses = num_subClasses + 1;
    subclassList{num_subClasses} = channels{iCh}.location.subclass;
    
end
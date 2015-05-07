function [sessionList] = getSessionsfromChannelDB( channels )
%
% usage: [sessionList] = getSessionsfromChannelDB( channels )
%
% function to extract the set of individual sessions from a channel DB
% structure
%
% INPUTS:
%   channels - a channel structure
%
% OUTPUTS:
%   sessionList - a cell array containing the list of sessions
%      included in the "channels" channel DB

sessionList{1} = channels{1}.session;
numSessions = 1;
for iCh = 1 : length(channels)
    
    if max(strcmp( channels{iCh}.session, sessionList ))
        % this session has already been counted
        continue;
    end
    
    numSessions = numSessions + 1;
    sessionList{numSessions} = channels{iCh}.session;
    
end
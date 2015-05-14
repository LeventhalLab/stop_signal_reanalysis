function newChannels = excludeChannels( channelParams, oldChannels )
%
% usage: newChannels = excludeChannels( channelParams, oldChannels )
%
% function to exclude channels with the given channelParams
%
% INPUTS:
%   channelParams - structure containing channel information to match
%       .task - 3 = stop-signal, 4 = go/nogo, -1 = either. Can be a vector
%           if, in the future more task options are available
%       .locationName - the brain region where the tetrode was.
%       .locationSubClass - location subclass
%       .subject - the implant ID# (of the form "IM-xxx")
%       .date - date in the format "yyyy-mm-dd'
%       .tetrode - recording site in the form 'Xnn', where "X" is "E", "R",
%           or "T", and nn is a number 01 - 18
%       .session - the session name of the form "Dxxyyyymmdd"; "D"
%           indicates one of Dan's sessions, xx is the original rat number, and
%          yyyymmdd is the date
%       .channelName - name of the individual channel in the form
%           "DxxyyyymmddXnn"; "D" indicates Dan's recording, xx is the
%           original rat number, yyyymmdd is the date, "Xnn" is the
%           tetrode/ref/EEG number
% for any numerical parameters (ie, task type, etc), a value of -1
% indicates that any value is acceptable. For strings, the string "any"
% means that any value is acceptable. Strings should be stored in cell
% arrays, so multiple values could be acceptable.
%
%   oldChannels - a channels structure
%
% OUTPUTS:
%   newChannels - a channel structure containing only the channels kept
%      after excluding any channels that fit channelParams

chList = extractChannels(channelParams, oldChannels);

channels_to_keep = ones(1, length(oldChannels));

channels_to_keep(chList) = 0;

newChannels = oldChannels(channels_to_keep == 1);
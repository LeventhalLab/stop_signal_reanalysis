function channelList = extractChannels( channelParams, channels, varargin )
%
% usage: channelList = extractChannels( channelParams, channels, varargin )
%
% function to extract channels with specific attributes (for example, by
% task, by brain region, number of a specific trial type, etc.)
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
%       .isValid - whether or not to accept only "valid" channels (
%           0 = accept invalid channels, 1 = accept only valid channels,
%           -1 = accept any channels)
%
% for any numerical parameters (ie, task type, etc), a value of -1
% indicates that any value is acceptable. For strings, the string "any"
% means that any value is acceptable. Strings should be stored in cell
% arrays, so multiple values could be acceptable.
%
% varargs: these have been included to account for any factors I didn't
% initially consider
%
% OUTPUTS:
%   a list of indices into the channels structure identifying channels that
%   match channelParams. Returns an empty array if no channels match the
%   requested parameters.

numChannels = length(channels);
numTargetChannels = 0;
channelList = [];

fieldList = fieldnames(channelParams);
% first, convert any string input into cells of strings
for iField = 1 : length(fieldList)
    
    if any(strcmpi(fieldList{iField}, {'task', 'isValid'}))
        continue;
    end
    if ~iscell(channelParams.(fieldList{iField}))
        channelParams.(fieldList{iField}) = ...
            cellstr(channelParams.(fieldList{iField}));
    end
    
end

for iChannel = 1 : numChannels
    
    % does the task match?
    if channelParams.task(1) > -1
        
        if ~max(channels{iChannel}.task == channelParams.task)
            % no match - don't extract this channel
            continue;
        end
        
    end
    
    % does the location match
    if ~strcmpi(channelParams.locationName{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.location.name, ...
                channelParams.locationName))
            continue;
        end
        
    end
    
    % does the location subclass match
    if ~strcmpi(channelParams.locationSubClass{1}, 'any')
        % change subclass to subClass for nico's data
        if ~max(strcmpi(channels{iChannel}.location.subclass, ...
                channelParams.locationSubClass))
            continue;
        end
        
    end
    
    if ~strcmpi(channelParams.subject{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.subject, ...
                channelParams.subject))
            continue;
        end
        
    end
    
    if ~strcmpi(channelParams.date{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.date, ...
                channelParams.date))
            continue;
        end
        
    end
    
    if ~strcmpi(channelParams.tetrode{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.tetrode.name, ...
                channelParams.tetrode))
            continue;
        end
        
    end
    
    if ~strcmpi(channelParams.session{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.session, ...
                channelParams.session))
            continue;
        end
        
    end
    
    if ~strcmpi(channelParams.channelName{1}, 'any')
        
        if ~max(strcmpi(channels{iChannel}.name, ...
                channelParams.channelName))
            continue;
        end
        
    end
    
    if channelParams.isValid > -1
        
        if isfield(channels{iChannel}, 'isValid')
            % perform this check only if the "isValid" field exists for
            % this channel structure; otherwise, accept all channels
            if channelParams.isValid ~= channels{iChannel}.isValid
                continue;
            end
        end
        
    end
    
    % this channel is a match
    numTargetChannels = numTargetChannels + 1;
    channelList(numTargetChannels) = iChannel;
    
end    % for iChannel = 1 : numChannels

if isempty(channelList)
    disp('No channels match these parameters');
end

            
        
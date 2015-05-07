function trIdx = extractTrials( t, trialEventParams )
%
% usage: trIdx = extractTrials( t, trialEventParams )
%
% function to find the indices of trials that meet the specifications of
% trialEventParams.
%
% INPUTS:
%   t - trial array
%   trialEventParams - a structure containing information necessary to extract
%      appropriate trials. If trialTypes is a cell array, trialEventParams is
%      an array of structures. the trialEventParams structure has the following
%      fields:
%           trialType - 
%              0 = any
%              1 = go
%              2 = no-go
%              3 = stop-signal
%              4 = stop converted to go
%              5 = any go trial
%              6 = selected at random (monte carlo)
%           correct - indicator of successful vs failed trial:
%              -1 = either, 0 = failure, 1 = correct
%           countsAsTrial - indicates whether or not to include only trials
%               that "counted": 0 - exclude "bad" trials, 1 - include "bad"
%               trials. NOTE: for GO trials, countsAsTrial excludes false
%               starts, wrong starts, but INCLUDES limited hold failures,
%               movement hold failures, and wrong target trials
%           tone - which tone played:
%               -1 - either tone; 1 - low tone; 2 - high tone
%           movementDirection: the direction the rat actually moved:
%               'Left', 'Right', or empty for either direction
%           outcome - string describing the trial's outcome:
%               'correct', 'wrongtarget', 'lhfail', 'failedstop',
%               'failednogo', 'any'
%
% OUTPUTS:
%   trIdx - vector containing the indices of trials that meet criteria

trialParamFields = fieldnames(trialEventParams);
numFields = length(trialParamFields);
numTrials = length(t);
tMarkers = false(numTrials, 1);

trList = true(numTrials, 1);

if isfield(t(1), 'trialType') && trialEventParams.trialType > 0
    
    trialType = [t.trialType]';
    
    if trialEventParams.trialType == 5 
        tMarkers = (trialType == 1 | trialType == 4);
    else
        tMarkers = (trialType == trialEventParams.trialType);
    end

else
    
    tMarkers(:,1) = true(numTrials, 1);
    
end

trList = trList & tMarkers;

for iField = 2 : numFields - 2

    if isfield(t(1), trialParamFields{iField})
        
        if trialEventParams.(trialParamFields{iField}) > -1
            
            for iTest = 1 : length(t)
                % this for loop is a work-around because occasionally some
                % values are missing for trial parameters - DL 1/15/2011
                if isempty(t(iTest).(trialParamFields{iField}))
                    t(iTest).(trialParamFields{iField}) = 0;
                end
            end
            
            fieldMarkers = [t.(trialParamFields{iField})]';
            tMarkers = (fieldMarkers == trialEventParams.(trialParamFields{iField}));
        else
            tMarkers = true(numTrials, 1);
        end   % if trialEventParams.(trialParamFields{iField}) > -1

        trList = trList & tMarkers;
        
    end
    
end    % end for iField...

trIdx = find(trList);
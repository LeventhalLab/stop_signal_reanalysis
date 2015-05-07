function [ trialEventParams ] = getTrialEventParams( trialTypes, varargin  )
% function to return a set of trial parameters for analysis given an input
% string describing the trial type
%
% INPUT:
%   trialTypes - a string or cell array of strings with descriptions of the
%       trials types to be analyzed (see possibilities in the "switch"
%       function below
%   eventString - name of events to be extracted - ie, "cueon", "nosein",
%       etc.
%
% varargins: used to specify subgroups of trials - for example, high or low
%   pitched tone, or moved left versus right
%       'tone': -1 (either), 1 (low), or 2 (high)
%       'movementdirection: -1 (either), 1 (left), or 2 (right)
%
% OUTPUT:
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

trialEventParams.trialType = 0;
trialEventParams.correct = 0;
trialEventParams.countsAsTrial = 0;
trialEventParams.tone = 0;
trialEventParams.movementDirection = 0;
trialEventParams.falseStart = 0;
trialEventParams.centerNP = 0;
trialEventParams.sideNP = 0;
trialEventParams.holdTooLong = 0;
trialEventParams.movementTooLong = 0;
trialEventParams.outcome = 0;
trialEventParams.event = '';
                 
switch lower(trialTypes)
    
    case 'montecarlo',
        trialEventParams.trialType = 6;
        trialEventParams.correct = [];
        trialEventParams.countsAsTrial = [];
        trialEventParams.tone = [];
        trialEventParams.movementDirection = [];
        trialEventParams.falseStart = [];
        trialEventParams.centerNP = [];
        trialEventParams.sideNP = [];
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'any';
        
    case {'allgo', 'all go'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = -1;
        trialEventParams.countsAsTrial = -1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = -1;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = -1;
        trialEventParams.movementTooLong = -1;
        trialEventParams.outcome = 'any';
        
    case {'completedgo', 'completed go'},   % this excludes LH and MH failures, at least for now
        trialEventParams.trialType = 5;
        trialEventParams.correct = -1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = 0;
        trialEventParams.movementTooLong = 0;
        trialEventParams.outcome = 'any';
        
    case {'correctgo', 'correct go'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = 1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = 0;
        trialEventParams.movementTooLong = 0;
        trialEventParams.outcome = 'correct';
        
    case {'correctgoright', 'correct go right'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = 1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = 2;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = 0;
        trialEventParams.movementTooLong = 0;
        trialEventParams.outcome = 'correct';


    case {'correctgoleft', 'correct go left'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = 1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = 1;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = 0;
        trialEventParams.movementTooLong = 0;
        trialEventParams.outcome = 'correct';


    case {'lhviol', 'lh viol', 'lhviolation', 'lh violation', 'lhfail', ...
            'lh fail'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = 0;
        trialEventParams.countsAsTrial = -1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = [];
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = [];
        trialEventParams.holdTooLong = 1;
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'lhfail';
        
    case {'wronggo', 'wrong go', 'wrong direction', 'wrongdirection', ...
            'wrong target', 'wrongtarget'},
        trialEventParams.trialType = 5;
        trialEventParams.correct = 0;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = 0;
        trialEventParams.movementTooLong = 0;
        trialEventParams.outcome = 'wrongtarget';
        
    case {'allstop', 'all stop', 'stop'},
        trialEventParams.trialType = 3;
        trialEventParams.correct = -1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = -1;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'any';
        
    case {'allnogo', 'all nogo', 'nogo'},
        trialEventParams.trialType = 2;
        trialEventParams.correct = -1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = -1;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'any';
        
    case {'correctstop', 'correct stop'},
        trialEventParams.trialType = 3;
        trialEventParams.correct = 1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = -1;
        trialEventParams.movementTooLong = -1;
        trialEventParams.outcome = 'correct';
        
    case {'correctnogo', 'correct nogo'},
        trialEventParams.trialType = 2;
        trialEventParams.correct = 1;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = [];
        trialEventParams.movementDirection = [];
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = [];
        trialEventParams.sideNP = [];
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'correct';
        
    case {'failedstop', 'failed stop'},
        trialEventParams.trialType = 3;
        trialEventParams.correct = 0;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = [];
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = [];
        trialEventParams.sideNP = [];
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'failedstop';
        
    case {'failednogo', 'failed nogo'},
        trialEventParams.trialType = 2;
        trialEventParams.correct = 0;
        trialEventParams.countsAsTrial = 1;
        trialEventParams.tone = [];
        trialEventParams.movementDirection = [];
        trialEventParams.falseStart = 0;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = [];
        trialEventParams.holdTooLong = [];
        trialEventParams.movementTooLong = [];
        trialEventParams.outcome = 'failedstop';
        
    case {'any', 'all'}
        trialEventParams.trialType = 0;
        trialEventParams.correct = -1;
        trialEventParams.countsAsTrial = -1;
        trialEventParams.tone = -1;
        trialEventParams.movementDirection = -1;
        trialEventParams.falseStart = -1;
        trialEventParams.centerNP = -1;
        trialEventParams.sideNP = -1;
        trialEventParams.holdTooLong = -1;
        trialEventParams.movementTooLong = -1;
        trialEventParams.outcome = 'any';
        
    otherwise,
        error('No match for trial type');
        
end

for iarg = 1 : 2 : nargin - 1
    
    switch lower(varargin{iarg})
        case 'movementdirection',
            trialEventParams.movementDirection = varargin{iarg + 1};
        case 'tone',
            trialEventParams.tone = varargin{iarg + 1};     
    end
    
end    % end for iarg...

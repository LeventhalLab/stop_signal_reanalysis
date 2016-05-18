function [RT, MT, sessionList] = collect_RT_MT_by_rat(channels, trialType)
%
% usage: [RT, MT, sessionList] = collect_RT_MT_by_rat(channels, trialType)
%
% INPUTS:
%   channels - a full channel structure
%   trialType - type of trials to use (see "getTrialEventParams" for a list
%       of possible trial types
%
% OUTPUTS:
%   RT - cell array containing a vector of RTs in trial order for the
%       selected trial type for all sessions
%   MT - cell array containing a vector of MTs in trial order for the
%       selected trial type for all sessions
%   sessionList - list of session names from the channel structure

trialEventParams = getTrialEventParams(trialType);
cp = initChanParams();
cp.locationSubClass = {'eeglam','ref','emg'};
cp.isValid = -1;
channels = excludeChannels(cp, channels);

cp = initChanParams();
cp.isValid = 1;
chList = extractChannels(cp, channels);
channels = channels(chList);

sessionList = getSessionsfromChannelDB(channels);
numSessions = length(sessionList);

RT = cell(1, numSessions);
MT = cell(1, numSessions);

for iSession = 1 : numSessions
%     iSession
    cp = initChanParams();
    cp.session = sessionList{iSession};
    chList = extractChannels(cp, channels);
    
    ch = channels(chList);
    c = ch{1};
    
    t = c.trials;
    t = extractTrials2(t, trialEventParams);
    
    numTrials = length(t);
    numValid_RT_Trials = 0;
    numValid_MT_Trials = 0;
    for iTrial = 1 : numTrials
        if isfield(t(iTrial).timing,'reactionTime')
            numValid_RT_Trials = numValid_RT_Trials + 1;
            RT{iSession}(numValid_RT_Trials) = ...
                t(iTrial).timing.reactionTime;
        end
        if isfield(t(iTrial).timing,'movementTime')
            numValid_MT_Trials = numValid_MT_Trials + 1;
            MT{iSession}(numValid_MT_Trials) = ...
                t(iTrial).timing.movementTime;
        end
    end
end
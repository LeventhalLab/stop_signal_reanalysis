function et_lfps = getEventTriggeredLFP(eventTimes, lfp, winStart, winDuration, eventDuration, varargin)
%
% usage: 
%
% INPUTS:
%   numEventTypes - vector of event timestamps
%   lfp - structure with the following features:
%       .lfp - the lfp
%       .Fs  - lfp sampling rate
%   winStart - vector containing timestamps at the start of each temporal
%       analysis window (in seconds)
%   winDuration - duration of each analysis window (in seconds)
%   eventDuration - amount of time to extract around each event; algorithm
%       assumes that the event is right in the middle (eventDuration/2
%       seconds before and eventDuration/2 seconds  after each event)
%
% VARARGINS:
%
% OUTPUTS:
%   et_lfps - event-triggered field potentials in an m x n matrix. Each
%       column is a vector of values from a single time window, individual
%       time windows are represented in rows

sampDuration = ceil(winDuration * lfp.Fs);
eventSamps   = ceil(eventDuration * lfp.Fs);
totalSamps   = length(lfp);

et_lfps = zeros(sampDuration, length(winStart));

% make sure lfp is a column vector
if length(lfp.lfp) > size(lfp.lfp, 1)
    lfp.lfp = lfp.lfp';
end
% make sure eventTimes is a column vector
if length(eventTimes) > size(eventTimes, 1)
    eventTimes = eventTimes';
end

% find events within all time windows of interest
validEvents = [];
for iWin = 1 : length(winStart)
    startSamp = round(winStart(iWin) * lfp.Fs);
    if startSamp < ceil(eventSamps / 2)
        continue;
    end
    endSamp   = startSamp + sampDuration;
    if endSamp > totalSamps - ceil(eventSamps / 2)
        continue;
    end
    validEvents = [validEvents; eventTimes(startSamp : endSamp)];
end

et_lfps = zeros(sampDuration, length(validEvents));

for iEvent = 1 : length(validEvents)
    
    eventStart = round(validEvents(iEvent) * lfp.Fs) - ceil(eventSamps / 2);
    eventEnd   = eventStart + eventSamps;
    et_lfps(:, iEvent) = lfp.lfp(eventStart : eventEnd);
    
end
function plvSurr = plv_surrogates( ch, ...
                                   sp_ts, ...
                                   trialType, ...
                                   eventList, ...
                                   eventWin, ...
                                   spikeWin, ...
                                   stepSize, ...
                                   surrogateSaveName, ...
                                   spikeLFP_surr_metadata, ...
                                   varargin )
% function 

hilbert_directory_1Hz = sprintf('/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins');
hilbert_directory_025Hz = sprintf('/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 025 Hz bins');


for iEventType = 1 : numEventTypes
    eventName = eventList{iEventType};
    event_ts = extractEvent_ts( validTrials, eventName );
    numEvents = length(event_ts);
    spikeLFP_metadata.numEvents = numEvents;
%         startSamps = round((event_ts + eventWin(1)) * metadata.Fs);
    eventStartTimes = event_ts + eventWin(1);
    
    for iStep = 1 : numSteps
%         tic
        validSpikes{iEventType, iStep} = [];
        winStartTimes = eventStartTimes + (iStep-1) * stepSize;
        winEndTimes   = winStartTimes + spikeWin;
        for i_singleEvent = 1 : numEvents
            validSpikes{iEventType, iStep} = [validSpikes{iEventType, iStep}; sp_ts(sp_ts > winStartTimes(i_singleEvent) & sp_ts < winEndTimes(i_singleEvent))];
        end
%         toc
        numSpikes(iEventType, iStep) = length(validSpikes{iEventType, iStep});

    end
        
end
function lfp_events = extractLFP(y, ch, trList, eventName, twin, Fs)
%
% INPUTS:
%   y - the full lfp
%   ch - single entry from the channel structure
%   trList - vector of trial structures
%   eventName - name of the trial event around which to extract LFPs

numTrials = length(trList);
samps_per_window = round(range(twin) * Fs);
lfp_events = zeros(numTrials, samps_per_window);

for iTr = 1 : numTrials
%     iTr
    if ~isfield(ch.trials(trList(iTr)).timestamps, eventName); continue; end
    event_ts = ch.trials(trList(iTr)).timestamps.(eventName);
    
    tlim = event_ts + twin;
    
    startSamp = round(tlim(1) * Fs);
    
    lfp_events(iTr, :) = y( startSamp : startSamp+samps_per_window-1 );
end
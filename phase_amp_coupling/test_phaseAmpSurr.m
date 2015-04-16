%%
numEvents = length(phaseAmp_metadata.eventList);
numPhaseFreq = length(phaseAmp_metadata.low_freqs);
numAmpFreq   = length(phaseAmp_metadata.high_freqs);
numSamps     = size(re_mrv, 4);
mrv_z = zeros(numEvents, numPhaseFreq, numAmpFreq, numSamps);
for iEvent = 1 : numEvents
    for i_f1 = 1 : numPhaseFreq
        for i_f2 = 1 : numAmpFreq
            for i_t = 1 : numSamps
                mrv_z(iEvent, i_f1, i_f2, i_t) = (abs(mrv(iEvent, i_f1, i_f2, i_t)) - surrogate_mean(i_f1, i_f2)) / surrogate_std(i_f1, i_f2);
                
            end
        end
    end
end

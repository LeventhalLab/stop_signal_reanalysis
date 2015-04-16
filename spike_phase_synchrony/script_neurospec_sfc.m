% script to call neurospec_sfc

eventList = {'cueOn','noseCenterIn','tone','noseCenterOut','noseSideIn','noseSideIn'};


ch = D22Ch_beta{745};
ts = D22spike{300}.timestamps.spike; 

trialType = 'correctgo';

[sfc, numSpikes] = calc_SFC_freqbands(ch, ts, trialType, eventList, [-1 1], 0.1, 0.05);
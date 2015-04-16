function calc_phaseRT_pvalue( ch, 
%
% function to bootstrap phase-RT correlations, per the analysis of Drewes
% and VanRullen, J Neurosci, 31(12):4698-708
%
% usage: 
%
% INPUTS: 
%   ch - a full channel DB structure for a single rat
%   

bitOrder = 'b';
trialType = 'correctgo';

hilbert_directory      = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/Hilbert transformed LFP 1 Hz bins';
phaseRTcorr_directory  = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT_correlations';
% RTcorr_plots_directory = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/phase_RT correlation plots';

implantID = implantID_from_ratID(ch.name(1:3));

subject_phaseRTcorr_directory = fullfile(phaseRTcorr_directory, [implantID '_phaseRTcorr']);
if ~exist(subject_phaseRTcorr_directory, 'dir')
    disp([subject_phaseRTcorr_directory ' not found.']);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% batchScript

% try
%     script_createRTphase_surrogates_b
% catch
% end
% 
try
    script_store_STOPsuccvfail_mrv_b;
catch
end

try
    script_store_STOPsuccvfail_mrv_surrogates;
catch
end

try
    script_plot_STOPsuccvfail_mrv_b
catch
end
% 
% try
%     script_windows_phaseAmp_coupling_2
% catch
% end
% 
% try
%     storeAnalyticSignals_025Hz_b
% catch
% end
% 
% try 
%     script_windows_phaseAmp_coupling_b;
% catch
% end
% 
% try
%     script_plot_all_trial_LFPs;
% catch
% end
% 
% try
%     script_phaseRThistograms;
% catch
% end
% 
% try
%     script_analyzeRTphasehists;
% catch
% end
% 
% try
%     storeAnalyticSignals
% catch
% end
% 
% try
%     script_Fries_sfc
% catch
% end
% 
% 
% % NEED TO LOOK AT DYNAMIC PHASE-AMPLITUDE COUPLING
% 
% % try
% %     script_storeMI
% % catch
% %     disp('error in script_storeMI')
% % end
% 
% try
%     script_plotMI
% catch
%     disp('error in script_plotMI')
% end
% 
% % try
% %     script_RTcorrelations
% % catch
% %     disp('error in script_RTcorrelations')
% % end
% 
% % CHECK TO SEE THAT THE CREATE_SURROGATES FUNCTION IS DOING WHAT I WANT...
% % try
% %     script_createRTphase_surrogates
% % catch
% %     disp('error in script_createRTphase_surrogates');
% % end
% 
% try
%     script_storeHilbertPowerComod
% catch
%     disp('error in script_storeHilbertPowerComod')
% end
% 
% 
% % analyses to consider:
% 
% %   CFC only during beta epochs?
% %   CFC on stop-success vs stop-failure, go-success vs go-failure
% %   CFC only during trials
% %   CFC immediately after the GO cue
% %   event-related coherence
% %   theta CFC at reward times?
% 
% %   expect dopamine to go up during initial approach, right after side-in?
% %       shift to theta as the modulating phase at nose-in?
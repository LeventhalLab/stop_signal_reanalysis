% batchScript

% try
%     storeAnalyticSignals
% catch
% end

% try
%     script_storeMI
% catch
%     disp('error in script_storeMI')
% end

% try
%     script_RTcorrelations
% catch
%     disp('error in script_RTcorrelations')
% end

try
    storeAnalyticSignals_Damien;
catch
    disp('error in storeAnalyticSignals_Damien');
end

try
    script_storeHilbertPowerComod
catch
    disp('error in script_storeHilbertPowerComod')
end
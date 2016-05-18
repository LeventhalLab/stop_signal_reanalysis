function new_t = extractTrials2(t, trialEventParams)

if trialEventParams.trialType > 0
    if trialEventParams.trialType == 5   % GO or STOP converted to GO
        numGOtrials = 0;
        for iTrial = 1 : length(t)
            if isfield(t(iTrial), 'trialType')
                if t(iTrial).trialType == 1 || t(iTrial).trialType == 4
                    numGOtrials = numGOtrials + 1;
                    t(numGOtrials) = t(iTrial);
                end
            else
                % this must be Greg's - only go trials exist
                numGOtrials = numGOtrials + 1;
                t(numGOtrials) = t(iTrial);
            end
        end    % end for iTrial...
            t = t(1 : numGOtrials);   % I went through this routine to make sure the order of the trials is preserved
    else
        if isfield(t(1),'trialType')
            t = t([t.trialType] == trialEventParams.trialType );
        elseif trialEventParams.trialType ~= 1
            % this must be one of Greg's trial structures, so no nogo or
            % stop trials. If go trials selected, "t" should include all
            % trials going forward
            new_t = [];
            return;
        end
    end
end    % end if trialEventParams.trialType > 0
        
trialParamFields = fieldnames(trialEventParams);
for iField = 2 : length(trialParamFields) - 2
    
    if isfield(t, trialParamFields{iField}) % some fields in DKL's data not present in GG's data
        if trialEventParams.(trialParamFields{iField}) > -1
            t  = t([t.(trialParamFields{iField})] == ...
                trialEventParams.(trialParamFields{iField})); 
        end % if trialEventParams.(trialParamFields{iField}) > -1
    end % if isfield(t, trialParamFields{iField})
end    % end for iField...

new_t = t;
function h_fig = phaseRTcorr_summary_sheet(chName, phaseRTcorr_metadata, RTphases, surrogate_metadata, mrl)

% adapted from VanRullen et al, Ongoing
% EEG phase..., Frontiers in Psychology, 2011; details in Drewes and
% Vanrullen, "This is the rhtyhm of your eyes: the phase of ongoing
% electroencephalogram oscillations modulates saccadic reaction time.", J
% Neurosci, 2011

% RTphases is a 3-dimensional (m x n x p) cell array. Each cell contains an
% array where each row contains the phases for a single trial around a
% specific event
%   m - number of events
%   n - number of center frequencies
%   p - RT quantile
%   
% mrl is an (m x n x p x q) array of bootstrapped mean resultant lengths
%   m - numFreqs
%   n - numEvents
%   p - each iteration of bootstrapping
%   q - time

trialType = phaseRTcorr_metadata.trialType;

%%
numEvents = length(phaseRTcorr_metadata.eventList);
numRTquantiles = length(phaseRTcorr_metadata.RTquantiles);
numFreqs  = length(phaseRTcorr_metadata.freqs);
num_t = size(RTphases{1,1,1}, 2);



numSurrogates = size(mrl, 3);

ITC = zeros(numEvents, numRTquantiles, num_t, numFreqs);   % inter-trial coherence (the "real" mrl)
surr_mean = zeros(numEvents, num_t, numFreqs);
p   = zeros(numEvents, numRTquantiles, num_t, numFreqs);   % p-value (where does ITC fall within the surrogate mrl distribution)

for iEvent = 1 : numEvents
    iEvent
    for iFreq = 1 : numFreqs
        for i_t = 1 : num_t
            surrogate_dist = squeeze(mrl(iFreq, iEvent, :, i_t));
            % get mean and standard deviations of surrogate distributions
            surr_mean(iEvent, i_t, iFreq) = mean(surrogate_dist);
            surr_std(iEvent, i_t, iFreq) = std(surrogate_dist);
            
            for iQuantile = 1 : numRTquantiles
        
            
                ph = RTphases{iEvent, iFreq, iQuantile}(:, i_t);
                r = sum(exp(1i*ph));
                ITC(iEvent, iQuantile, i_t, iFreq) = abs(r) / length(ph);
                
                p(iEvent, iQuantile, i_t, iFreq) = normpdf(ITC(iEvent, iQuantile, i_t, iFreq), ...
                                                           surr_mean(iEvent, i_t, iFreq), ...
                                                           surr_std(iEvent, i_t, iFreq));
                
            end
        end
    end
end

%%
% units are centimeters
figProps.ltMargin = 1.5;
figProps.rtMargin = 1.5;
figProps.topMargin = 2;
figProps.botMargin = 1;

figProps.width = 11 * 2.54;
figProps.height = 8.5 * 2.54;

figProps.n = numRTquantiles;
figProps.m = numEvents;

figProps.colSpacing = ones(1, figProps.n-1) * 0.1;
figProps.rowSpacing = ones(1, figProps.m-1) * 0.15;

totFigWidth = figProps.width - figProps.ltMargin - figProps.rtMargin;
totFigHeight = figProps.height - figProps.topMargin - figProps.botMargin;
figProps.panelWidth  = ones(1, figProps.n) * (totFigWidth - sum(figProps.colSpacing)) / figProps.n;
figProps.panelHeight = ones(1, figProps.m) * (totFigHeight - sum(figProps.rowSpacing)) / figProps.m;

h_fig = zeros(1,3);
[h_fig(1), h_axes] = createFigPanels4(figProps);
h_figAxes = createFigAxes(h_fig);
axes(h_figAxes);
textString{1} = ['mrl by RT quantile, ', chName, ', ' trialType];
text(0.05, 0.95, textString);

collim = [0 1];
t = linspace(-1,1,num_t);
for iEvent = 1 : numEvents
    for iQuantile = 1 : numRTquantiles
        
        axes(h_axes(iQuantile, iEvent));
        toplot = squeeze(ITC(iEvent, iQuantile, :, :))';
        
        imagesc(t, phaseRTcorr_metadata.freqs, toplot);
        
        set(gca,'ydir','normal', 'clim', collim);
        if iQuantile == 1
            title(phaseRTcorr_metadata.eventList{iEvent});
        end
        if iEvent > 1
            set(gca,'yticklabel',{})
        else
            ylabel(num2str(phaseRTcorr_metadata.RTquantiles(iQuantile)));
        end
        if iQuantile < numRTquantiles
            set(gca,'xticklabel',{})
        end
        
    end
    
end

[h_fig(2), h_axes] = createFigPanels4(figProps);
h_figAxes = createFigAxes(h_fig);
axes(h_figAxes);
for iEvent = 1 : numEvents
    for iQuantile = 1 : numRTquantiles
        
        axes(h_axes(iQuantile, iEvent));
        toplot = squeeze(p(iEvent, iQuantile, :, :))';
        
        imagesc(t, phaseRTcorr_metadata.freqs, toplot);
        
        set(gca,'ydir','normal', 'clim', collim);
        if iQuantile == 1
            title(phaseRTcorr_metadata.eventList{iEvent});
        end
        if iEvent > 1
            set(gca,'yticklabel',{})
        else
            ylabel(num2str(phaseRTcorr_metadata.RTquantiles(iQuantile)));
        end
        if iQuantile < numRTquantiles
            set(gca,'xticklabel',{})
        end
        
    end
end

figProps.m = 1;
figProps.botMargin = 3; figProps.topMargin = 4;
figProps.panelHeight = ones(1, figProps.m) * (totFigHeight - sum(figProps.rowSpacing)) / figProps.m;
[h_fig(3), h_axes] = createFigPanels4(figProps);
h_figAxes = createFigAxes(h_fig);
axes(h_figAxes);
for iEvent = 1 : numEvents
        
        axes(h_axes(iQuantile, iEvent));
        toplot = squeeze(surr_mean(iEvent, :, :))';
        
        imagesc(t, phaseRTcorr_metadata.freqs, toplot);
        
        set(gca,'ydir','normal', 'clim', collim);
        if iQuantile == 1
            title(phaseRTcorr_metadata.eventList{iEvent});
        end
        if iEvent > 1
            set(gca,'yticklabel',{})
        else
            ylabel(num2str(phaseRTcorr_metadata.RTquantiles(iQuantile)));
        end
        if iQuantile < numRTquantiles
            set(gca,'xticklabel',{})
        end
        
end
% 
% 
%         surr_to_plot = squeeze(surr_mean(iEvent, :, :))';
%         figure
%         imagesc(t, phaseRTcorr_metadata.freqs, surr_to_plot);
%         colorbar
%         set(gca,'ydir','normal','clim',collim);
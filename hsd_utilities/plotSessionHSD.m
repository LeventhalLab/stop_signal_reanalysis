function plotSessionHSD( varargin )
% function to plot every tetrode, reference, and EEG wire so that the good
% wires can be assessed for a given session

ylim = [-10 10];
pn = fullfile('/Volumes','dan-nas1-1','Recording data');
Fs = 31250;
timeLimits = [100 120];
fn = [];

if nargin / 2 ~= floor(nargin / 2)
    % not an even number of arguments
    disp('invalid number of input arguments')
    return
end    % end if nargin...

for iarg = 1 : 2 : nargin    % assumes an even number of varargs
    
    switch lower(varargin{iarg})
        
        case 'filename',
            fn = varargin{iarg + 1};
            fullName = fn;
            
        case 'fs',
            Fs = varargin{iarg + 1};
            
        case 'timelimits',
            timeLimits = varargin{iarg + 1};
           
        case 'ylim',
            ybound = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...

if isempty(fn)
    cd(pn)

    [fn, pn, FilterIndex] = uigetfile({'*.hsd'; '*.hsdw'});

    fullName = [pn filesep fn];
    
end    % end if isempty(fn)

% make sure that the time window to plot fits within the duration of the
% actual .hsd (or .hsdf) file
hsdDuration = getHSDlength( fullName );
sampleDuration = range(timeLimits);
if timeLimits(2) > hsdDuration
    timeLimits(2) = hsdDuration;
    if sampleDuration > hsdDuration
        timeLimits(1) = 0;
    else
        timeLimits(1) = hsdDuration - sampleDuration;
    end
end
% above code makes sure that timeLimits(2) is less than the total length in
% seconds of the .hsd (.hsdf) file. Then, it preserves the width of the
% time window to plot by going backwards from the end of the file (that is,
% if 2 seconds were to be plotted, it still plots 2 seconds by going back
% from the end of the file). If the width of the window to be plotted was
% greater than the entire length of the file, the start time is just the
% beginning of the file. -DL 4/12/2010
disp(['Start time: ' num2str(timeLimits(1)) ', End time: ' num2str(timeLimits(2))]);

for iTetrode = 1 : 18
    
    fhandle = figure;
    set(gcf, 'Name', ['tetrode ' num2str(iTetrode)]);
    
    wireNums = getTetrodeWires( iTetrode );
    
    plotHSD('filename',fullName,'plotwires',wireNums,'ylim',ylim,...
        'fs',Fs,'timelimits',timeLimits,'fhandle',fhandle);
        
        
    
end    % end for iTetrode

clear wireNums

for iRef = 1 : 3
    
    fhandle = figure;
    set(gcf, 'Name', ['reference ' num2str(iRef)]);
    
    wireNums = getRefWires( iRef );
    
    plotHSD('filename',fullName,'plotwires',wireNums,'ylim',ylim,...
        'fs',Fs,'timelimits',timeLimits,'fhandle',fhandle);
    
end    % end for iRef

for iEEG = 1 : 3
    
    fhandle = figure;
    set(gcf,'Name', ['EEG ' num2str(iEEG)]);
    
    wireNum = getEEGWire( iEEG );
    
    plotHSD('filename',fullName,'plotwires',wireNum,'ylim',ylim,...
        'fs',Fs,'timelimits',timeLimits,'fhandle',fhandle);
    
end    % end for iEEG
    
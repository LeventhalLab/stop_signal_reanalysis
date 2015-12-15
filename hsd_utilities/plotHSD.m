function [plotData] = plotHSD( varargin )
%
% function to load and plot a tetrode's worth of raw/filtered data
% plotWires is the numbers of the channels to be plotted.
%
% IM 164 hsdw files
% pn = '/Volumes/dan/Recording data/IM-164_DL-22/IM-164_DL-22_hsd/IM164_HSDw/IM164_HSDw_20091113/';
%
% function arguments:
% 'filename' - name of .hsd or .hsdw file to read in
% 'numwires' - the number of wires stored in the file
% 'fs' = sampling frequency
% 'timelimits' - 1 x 2 array with the start and end times to plot
% 'plotwires' - 1 x n array with the indices of the wires to plot
% 'ylim' - 1 x 2 array with the y axis min and max (in mV)
%
% OUTPUT:
%   plotData is a m x n array, where m is the number of plots made, and n
%   is the length of the plots (in samples).

if nargin / 2 ~= floor(nargin / 2)
    % not an even number of arguments
    disp('invalid number of input arguments')
    return
end    % end if nargin...

numWires = 4;
Fs = 31250;
timeLimits = [16 18];
plotWires = [1 : 4];
plotWires2 = 0;
%defaultPath = '/Volumes/dan/Recording data/IM-164_DL-22/IM-164_DL-22 screw turn hsds/';
%defaultPath = '/Volumes/dan/Recording data/IM-166_DL-20/IM-166_DL-20 screw turn hsds/';
%defaultPath = '/Volumes/dan/Recording data/IM-164_DL-22/IM-164_DL-22_hsd';
defaultPath = cd;
fhandle = [];
axis_h = [];
convertToVolts = 1;
isVisible = 'on';

ybound = [];
fn = '';


for iarg = 1 : 2 : nargin    % assumes an even number of varargs
    
    switch lower(varargin{iarg})
        
        case 'filename',
            fn = varargin{iarg + 1};
            fullName = fn;
        case 'numwires',
            numWires = varargin{iarg + 1};
            plotWires = [1 : numWires];
        case 'fs',
            Fs = varargin{iarg + 1};
        case 'timelimits',
            timeLimits = varargin{iarg + 1};
        case 'plotwires',
            plotWires2 = varargin{iarg + 1};
        case 'ylim',
            ybound = varargin{iarg + 1}; 
        case 'fhandle',
            fhandle = varargin{iarg + 1};
        case 'axeshandle',
            axis_h = varargin{iarg + 1};
        case 'converttovolts',
            convertToVolts = varargin{iarg + 1};
        case 'isvisible',
            isVisible = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...

if plotWires2
    plotWires = plotWires2;
    numWires = length(plotWires);
end    % end if plotWires2

if ~isempty(axis_h)
    if length(axis_h) ~= numWires
        disp('Mismatch between number of wires to plot and number of axes handles passed into plotHSD.');
        return;
    end
end

if isempty(fn)
    cd(defaultPath)

    [fn, pn, FilterIndex] = uigetfile({'*.hsd'; '*.hsdw'});

    fullName = [pn filesep fn];
    
end    % end if isempty(fn)

% check if the .hsd(w) file has a header
hasHeader = checkHSDheader( fullName );

if hasHeader
    % read in appropriate metadata
    hsdHeader = getHSDHeader( fullName );
    dataOffset = hsdHeader.dataOffset;
    numChannels = hsdHeader.main.num_channels;   % total number of wires (channels) in the hsd file
    Fs = hsdHeader.main.sampling_rate;
else
    % there was no header; data offset is zero
    dataOffset = 0;
end

rawData = readHSD( fullName, numChannels, dataOffset, Fs, timeLimits );

if convertToVolts
    rawData = int2volt( rawData );
end

if isempty(fhandle) && isempty(axis_h)
    figure;
    set(gcf, 'Name', [fn ' wires ' num2str(plotWires)]);
else
    if isempty(axis_h)
        figure(fhandle);
    end
end

x = 1/Fs : 1 / Fs : range(timeLimits);
plotData = rawData(plotWires, :);

for iWire = 1 : length(plotWires)
    if isempty(axis_h)
        subplot(length(plotWires), 1, iWire);
    else
        axes(axis_h(iWire));%, 'visible', isVisible);
    end

    plotRow = plotWires(iWire);
    
    plot(x, rawData(plotRow,:));
    text('units','normalized','position',[0.05 0.8],'string',num2str(plotRow));
    if ~isempty(ybound)
        set(gca, 'ylim', ybound);
    end

end

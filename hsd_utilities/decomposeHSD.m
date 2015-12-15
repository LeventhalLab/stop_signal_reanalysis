function [filenames] = decomposeHSD( varargin )

% function to split a .hsd file into its component wires, and save the
% files into the same folder from which the .hsd file was read

% varargs:
%   'filename' - name of .hsd file to decompose
%   'outpath' - path to output the new files to (default the same folder as
%       input)
%   'numwires' - the number of individual recording channels
%   'wirestosave' - array with the numbers of the wires to save; default is
%           all of them
%   'datatype' - format of data; 'int16' is default

% OUTPUTS:
%   filenames - 


fn = '';
numWires = 81;
wiresToSave = [1 : 81];
dataType = 'int16';
samples_per_block = 200000;
outPath = '';

%defaultPath = fullfile('/Volumes', 'dan', 'Recording data', ...
%    'IM-164_DL-22', 'IM-164_DL-22 screw turn hsds');
%cd(defaultPath);

for iarg = 1 : 2 : nargin
    
    switch lower(varargin{iarg})
        
        case 'filename',
            fn = varargin{iarg + 1};
            [pn, fn, ext, versn] = fileparts( fn );
            fn_ext = [fn ext];
            if ~isempty(pn)
                cd(pn);
            end
        case 'outpath',
            outPath = varargin{iarg + 1};
        case 'numwires',
            numWires = varargin{iarg + 1};
        case 'wirestosave',
            wiresToSave = varargin{iarg + 1};
        case 'datatype',
            dataType = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...

if isempty(fn)    % a filename wasn't specified
    
    
    [fn, pn] = uigetfile('*.hsd');
    
    fullname = [pn filesep fn];
    
	[pn, fn, ext, versn] = fileparts( fullname );
    fn_ext = [fn ext];
    cd( pn );
    
end    % end if isempty(fn)



hasHeader = checkHSDheader( fn_ext );

if hasHeader
    % read in metadata
    hsdHeader = getHSDHeader( fn_ext );
    dataOffset = hsdHeader.dataOffset;
    numWires = hsdHeader.main.num_channels;
else
    % there was no header; data offset is zero
    dataOffset = 0;
end

fileInfo = dir( fn_ext );
dataBytes = fileInfo.bytes - dataOffset;
bytes_per_wire = dataBytes / numWires;
bytes_per_sample = getBytesPerSample( dataType );
samples_per_wire = bytes_per_wire / bytes_per_sample;   % the number of samples total on a single wire

if isempty(outPath)
    outPath = pn;
end    % if a separate path has not been specified for the output files, set outPath to the .hsd path

% create array of output filenames
for iWire = 1 : length(wiresToSave)
    
    if wiresToSave(iWire) < 10
        outputFn(iWire,:) = [fn '_w0' num2str(wiresToSave(iWire)) '.wire.hsd'];
    else
        outputFn(iWire,:) = [fn '_w' num2str(wiresToSave(iWire)) '.wire.hsd'];
    end

    fout(iWire) = fopen(fullfile(outPath, outputFn(iWire,:)), 'w');
    
end    % end for iWire...
            
numBlocks = ceil(samples_per_wire / samples_per_block);

fin = fopen(fn_ext, 'r');
if fin == -1
    keyboard;
end
fseek(fin, dataOffset * bytes_per_sample, 'bof');

rawData = zeros(numWires, samples_per_block);
rawData = eval([dataType '(rawData);']);


for iBlock = 1 : numBlocks
%    disp(sprintf('%d of %d blocks...', iBlock, numBlocks));
    [rawData, num_elements_read] = ...
        fread(fin, [numWires, samples_per_block], dataType, 0, 'b');
    
    for iWire = 1 : length(wiresToSave)
        try
        	fwrite(fout(iWire), rawData(wiresToSave(iWire), :),dataType);
		catch
			keyboard;
		end
        
    end   % end for iWire...
end

fclose(fin);
filenames = cellstr(outputFn);

for iWire = 1 : length(wiresToSave)
    
    fclose(fout(iWire));
    
end    % end for iWire...


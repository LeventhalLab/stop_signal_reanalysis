function [maxTime, varargout] = getHSDlength( varargin )
%
% usage: [maxTime] = getHSDlength( varargin )
%
% function to get the length of an HSD file in seconds. If a filename is
% not provided, a get file dialog box is opened
%
% varargin:
%   filename - a string or cell array of strings containing file names
%
% OUTPUT:
%   maxTime - the length of the recording contained in the given file in
%       seconds. If multiple filenames are provided/selected, maxTime is a
%       vector containing the length in seconds of each file

maxTime = [];
dataType = 'int16';

if (nargin / 2) ~= round(nargin / 2)
    disp('function "getHSDlength" requires an even number of input arguments');
    return;
end

fn = [];
numChannels = 0;
for iarg = 1 : 2 : nargin
    
    switch lower(varargin{iarg})
        
        case 'filename',
            fn = varargin{iarg + 1};
            [pn, fn, ext] = fileparts( fn );
            if isempty(pn)
                pn = cd;
            end
            fn = [fn ext];
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'fs',
            Fs = varargin{iarg + 1};
            % if sampling rate is recorded in the header, the header Fs
            % will be used. If there is no header, Fs needs to be supplied
            % as a variable input argument.
        case 'numchannels',
            numChannels = varargin{iarg + 1};
            
    end    % end switch
    
end     % end for iarg...

bytes_per_sample = getBytesPerSample( dataType );

if isempty(fn)
    
    [fn, pn] = uigetfile({'*.hsd;*.hsdf', '.hsd or lfp files'}, ...
        'multiselect', 'on');
end

if ~iscell(fn)
    fn = cellstr(fn);
end

numFiles = length(fn);

maxTime = zeros(1, numFiles);
duration_in_samples = zeros(1, numFiles);

for iFile = 1 : numFiles

    fn{iFile} = fullfile(pn, fn{iFile});
    % check to make sure that fn{iFile} exists
    if ~exist(fn{iFile}, 'file')
        disp(['error in function "getHSDlength": file ' fn{iFile} ' does not exist.'])
        maxTime(iFile) = 0;
        continue;
    end
    
    f_info = dir( fn{iFile} );
    
    if checkHSDheader(fn{iFile})
        hsd_header = getHSDHeader( fn{iFile} );

        if strcmp(fn{iFile}(length(fn{iFile})), 'f')
            % this is an LFP (.hsdf) file
            if isfield(hsd_header.main, 'downsampled_rate')
                Fs = hsd_header.main.downsampled_rate;
                if Fs < 100
                    Fs = lfpFs( hsd_header );
                end
                % enough workarounds to get the sampling rate?
            else
                Fs = lfpFs(hsd_header);
            end
        else
            Fs = hsd_header.main.sampling_rate;
        end
        dataOffset = hsd_header.dataOffset;
        numChannels = hsd_header.main.num_channels;
        
    else
        % no header - assume no data offset
        dataOffset = 0;
        % if numChannels not defined as a varargin, assume a single channel
        if numChannels == 0; numChannels = 1; end
    end
    
    if ~exist('Fs','var')
        disp(['sampling rate not found for ' fn{iFile}]);
        continue;
    end
    
    numBytes = f_info.bytes;
    numDataBytes = numBytes - dataOffset;
    bytes_per_channel = numDataBytes / numChannels;
    samples_per_channel = bytes_per_channel / bytes_per_sample;
    
    maxTime( iFile ) = samples_per_channel / Fs;
    duration_in_samples(iFile) = samples_per_channel;
end    % end for iFile

varargout(1) = {duration_in_samples};
function data_out = volt2int( data_in, varargin )

% function to convert data from voltages into int16 (or other desired data
% type)

% input arguments:
% data_in - voltage to be converted into an integer, given in mV
%   datatype - data type to be converted into (default int16)
%   range - 1 x 2 array with lower and upper voltage range on daq card
%   (default [-10 10])
%   gain - amplitude gain on amplifiers (default 2000)

% output arguments:
%   data_out - signal voltages in millivolts

% program could be touched up to verify that variable input arguments are
% of types compatible with the m-file (ie, that range is a 1 x 2 numerical
% array)

dataType = 'int16';
voltRange = [-10 10];
ampGain = 2000;

for iarg = 1 : 2 : nargin - 1
    
    switch lower(varargin{iarg})
        
        case 'datatype',
            dataType = varargin{iarg + 1};
            
        case 'range',
            voltRange = varargin{iarg + 1};
            
        case 'gain',
            ampGain = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...
numBits = 8 * getBytesPerSample(dataType);

% above is to account for positive/negative numbers

intScale = ampGain / (range(voltRange) / 2^numBits);

data_out = eval([dataType '(data_in * intScale / 1000);']);
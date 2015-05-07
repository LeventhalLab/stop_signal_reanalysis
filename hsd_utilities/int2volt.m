function data_out = int2volt( data_in, varargin )
%
% usage: data_out = int2volt( data_in, varargin )
%
% function to convert data from int16 into voltages (output in mV)
% input arguments:
% data_in - integer data to be converted into voltages 
%
% varargs:
%   'precision' - number of bits of precision. default 16 (int16)
%   'range' - min and max voltages on daq card in V (default [-10 10])
%   'gain' - amplifier gain - default 2000
%
% output arguments:
%   data_out - signal voltages in millivolts
%
% program could be touched up to verify that variable input arguments are
% of types compatible with the m-file (ie, that range is a 1 x 2 numerical
% array)

numBits = 16;
voltRange = [-10 10];
ampGain = 2000;

for iarg = 1 : 2 : nargin - 1
    
    switch lower(varargin{iarg})
        
        case 'precision',
            numBits = varargin{iarg + 1};
        case 'range',
            voltRange = varargin{iarg + 1}; 
        case 'gain',
            ampGain = varargin{iarg + 1};
            
    end    % end switch
    
end    % end for iarg...

voltScale = (range(voltRange) / 2^numBits) / ampGain;

data_out = data_in * voltScale * 1000;
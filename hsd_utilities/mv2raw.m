function [ rawValue ] = mv2raw( mv, varargin )
%mv2raw
%   function to convert a value in millivolts to a raw data acquisition
%   integer value.

% INPUT:
%   mv - voltage in mv

% varargins:
%   'datatype' - the original recording data type; default 'int16'
%   'ampgain' - the amplifier gain from the recording; default 2000
%   'voltrange' - 2 element vector with the min and max voltages recorded
%       by the daq card. Default [-10 10]
% OUTPUT:
%   rawValue - the "raw data" value corresponding to the voltage in mv fed
%   to the function

dataType = 'int16';
ampGain = 2000;
voltRange = [-10 10];

for iarg = 1 : 2 : nargin - 1
    
    switch lower(varargin{iarg})
        
        case 'datatype',
            dataType = varargin{iarg + 1};
        case 'ampgain',
            ampGain = varargin{iarg + 1};
        case 'voltrange',
            voltRange = varargin{iarg + 1};
            
    end
    
end

numBits = 8 * getBytesPerSample( dataType );
voltScale = (range(voltRange) / 2^numBits) / ampGain;

rawValue = mv / voltScale / 1000;

end
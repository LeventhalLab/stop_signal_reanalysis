function numBytes = getBytesPerSample( datatype )

% function to return the number of bytes in a given data type
% if datatype is not recognized, the function returns 0

switch datatype
    
    case {'int8', 'uint8'},
        numBytes = 1;
    
    case {'int16', 'uint16'},
        numBytes = 2;
        
    case {'int32', 'uint32', 'single'},
        numBytes = 4;
        
    case {'int64', 'uint64', 'double'},
        numBytes = 8;
        
    otherwise,
        numBytes = 0;
        
end    % end switch
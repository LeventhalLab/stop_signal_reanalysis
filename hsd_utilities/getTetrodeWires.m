function [ wireNums ] = getTetrodeWires( iTetrode )

% function returns the wire numbers for a given tetrode in the Berke lab 
% standard  21 - tetrode drive

if iTetrode < 7
    
    wireNums = [(iTetrode - 1) * 4 + 1 : iTetrode * 4];
    
elseif iTetrode < 13
    
    wireNums = [iTetrode * 4 : iTetrode * 4 + 3];
    
else
    
    wireNums = [iTetrode * 4 + 3 : iTetrode * 4 + 6];
        
end    % end if iTetrode < 7
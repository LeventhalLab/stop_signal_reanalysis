function [ wireNum ] = getEEGWire( iEEG )

% function returns the wire number for a given EEG in the Berke lab 
% standard  21 - tetrode drive

switch iEEG
    
    case 1,
        
        wireNum = 27;
        
    case 2,
        
        wireNum = 54;
        
    case 3,
        
        wireNum = 81;
        
end    % end switch
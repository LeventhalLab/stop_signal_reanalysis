function [ wireNums ] = getRefWires( iRef )

% function returns the wire numbers for a given ref (stereotrode) in the 
% Berke lab standard  21 - tetrode drive

switch iRef
    
    case 1,
        
        wireNums = [25 : 26];
    
    case 2,
    
        wireNums = [52 : 53];
    
    case 3,
    
        wireNums = [79 : 80];
        
end    % end switch
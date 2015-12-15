function storeHeaders( hsdDir, targetDir )
%
% function to write the headers of a group of .hsd files into .hsdh files
% for easy storage of channel information that has been updated
%
% usage: storeHeaders( hsdDir, targetDir )
%
% INPUTS:
%   hsdDir - directory in which to find .hsd files
%   targetDir - directory in which to place .hsdh files
%

cd(hsdDir);
hsdList = dir('*.hsd');
numHSD = length(hsdList);

for iHSD = 1 : numHSD
    hsdName = hsdList(iHSD).name;
    
    header = getHSDHeader(hsdName);
    
    hsdh_name = fullfile(targetDir, [hsdName 'h']);
    
    writeHSDheader(hsdh_name, header);
    
end
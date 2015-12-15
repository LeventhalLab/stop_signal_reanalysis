function create_hsdh_from_folder(startDir, targetDir)
%
% usage: create_hsdh_from_folder(startDir, targetDir)
%
% function to take all the .hsd files contained within a given directory,
% extract their headers, and write just those headers into separate .hsdh
% files so that the metadata from a recording session can be moved around
% easilty without dragging gigabytes of data around as well.
%
% INPUTS:
%   startDir - folder in which the .hsd files are contained from which to
%       extract the headers
%   targetDir - folder in which to store the .hsdh files
%

cd(startDir)

fileList = dir('*.hsd');

for iFile = 1 : length(fileList)
    fnames{iFile} = fullfile(startDir, fileList(iFile).name);
end

make_hsdh_from_hsd(fnames, targetDir);
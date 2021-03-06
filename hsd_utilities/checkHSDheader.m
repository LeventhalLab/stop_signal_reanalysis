function [hasHeader] = checkHSDheader( filename )
%
% usage: [hasHeader] = checkHSDheader( filename )
%
% function to check if a .hsd file has a header or not. checks by seeing if
% the first four characters in the file are "HSD ". If not, then there is
% no header

fid = fopen(filename, 'r');

testString = 'HSD ';
readString = fread(fid, 4, '*char');

if strcmp(readString', testString)
    hasHeader = 1;
else
    hasHeader = 0;
end

fclose(fid);
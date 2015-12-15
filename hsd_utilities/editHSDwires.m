function editHSDwires( fn, editFeature, wireNums, newValues, varargin )
% usage: editHSDwires( fn, editFeature, wireNums, newValues, varargin )
% function to revise the header of a .hsd file to change wire information
%
% INPUTS:
%   fn - filename of the .hsd file; if empty, will prompt for a file
%   editFeature - feature that is to be changed. Options are:
%       'original_number'
%       'good'
%       'channel_type'
%       'channel_number'
%       'wire_number'
%       'bank'
%       'gain'
%       'low_cut'
%       'high_cut'
%       'name' - format 't01_1', 'r01_1'
%   wireNums - vector of wire numbers (1-81) for a 21-tetrode drive on
%       which changes will be made
%   newValues - vector of new values; must have the same length as wireNums
%       If editFeature is 'name', newValues should be a cell array of
%       strings
%
%   varargins:
%       'newfn' - new file to write to

editFeature = lower(editFeature);

if length(wireNums) ~= length(newValues)
    disp('wire number and new value vectors must have the same length');
    return;
end

if isempty( fn )
    
    % add code to pull up a get file gui
    
end

write_fn = fn;
[pathstr, fname, f_ext, versn] = fileparts( fn );

for iarg = 1 : 2 : nargin - 4
    
    switch(lower(varargin{iarg}))
        case 'newfn',
            write_fn = varargin{iarg + 1};
            if isempty(findstr(write_fn,'.'))
                % no file extension was included; make the file extension
                % the same as for the original file
            	write_fn = [write_fn f_ext];
            end
            
    end    % end switch
    
end    % end for iarg...


hsd_header = getHSDHeader( fn );

for i = 1 : length(newValues)
    
    if strcmp(editFeature, 'name')
        space_pad(1 : (64 - length(newValues{i}))) = ' ';
        newName = [newValues{i} space_pad];
        hsd_header.channel(wireNums(i)).name = newName;
    else
        evalStr = ['hsd_header.channel(wireNums(i)).' editFeature ' = ' num2str(newValues(i)) ';'];
        eval(evalStr);
    end    % end if strcmp(editFeature, 'name')
    
end    % end for i...

writeHSDheader(write_fn, hsd_header);


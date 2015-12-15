function make_hsdh_from_hsd( hsdFile, targetDir )
%
% usage: make_hsdh_from_hsd( hsdFile, targetDir )
%
% INPUTS:
%   hsdFile - name of the .hsd file(s) from which to pull the header(s)
%   targetDir - directory(s) in which to write the header(s) as a .hsdh
%       file(s)

if ~iscell(hsdFile)
    hsdFile = cellstr(hsdFile);    
end
if ~iscell(targetDir)
    targetDir = cellstr(targetDir);
end

if length(targetDir) ~= length(hsdFile)
    temp = targetDir{1};
    clear targetDir
    for i = 1 : length(hsdFile)
        targetDir{i} = temp;
    end
end
for ihsd = 1 : length(hsdFile)
    hsdHeader = getHSDHeader( hsdFile{ihsd} );

    [~, fname, ext, ~] = fileparts( hsdFile{ihsd} );
    
    hsdh_fn = fullfile(targetDir{ihsd}, [fname ext 'h']);

    writeHSDheader( hsdh_fn, hsdHeader);
end
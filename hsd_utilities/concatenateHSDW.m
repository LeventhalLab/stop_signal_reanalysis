function concatenateHSDW(filenames, outFile)
%
% usage: concatenateHSDW(filenames, outFile)
%
% INPUTS:
%   filenames - cell array of filenames

blockSize = 100000;

if exist(outFile,'file')
    overWriteOld = questdlg([outFile ' already exists. Overwrite?']);
    switch overWriteOld
        case 'No',
            [fn, pn. filtIndex] = uiputfile('*.hsdw');
            outFile = fullfile(pn, fn);
        case 'Cancel',
            return;
    end
end
            
    
fout = fopen(outFile, 'w');
for iFile = 1 : length(filenames)
    
    fin = fopen(filenames{iFile});
    if fin == -1
        disp(['error opening ' filenames{iFile}]);
        fclose(fout);
        return;
    end
    while ~feof(fin)
        hsdw = fread(fin, blockSize, 'int16',0,'b');
        fwrite(fout, hsdw,'int16',0,'b');
    end
    fclose(fin);
end

fclose(fout);
function concatenateHSD(filenames, outFile)

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
            
hsd_header = getHSDHeader(filenames{1});
numChannels = hsd_header.main.num_channels;
dataOffset = hsd_header.dataOffset;

writeHSDheader(outFile, hsd_header);

fout = fopen(outFile, 'r+');
fseek(fout, dataOffset, 'bof');
for iFile = 1 : length(filenames)
    
    fin = fopen(filenames{iFile});
    if fin == -1
        disp(['error opening ' filenames{iFile}]);
        fclose(fout);
        return;
    end
    curBlock = 0;
    fseek(fin, dataOffset, 'bof');
    while ~feof(fin)
        curBlock = curBlock + 1
        hsd = fread(fin, [numChannels, blockSize], 'int16',0,'b');
        fwrite(fout, hsd,'int16',0,'b');
    end
    fclose(fin);
end

fclose(fout);
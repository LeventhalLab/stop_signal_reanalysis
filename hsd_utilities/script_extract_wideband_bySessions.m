% script to automate extracting waveforms across multiple sessions

sessionsRoot = fullfile('/Volumes','dan-nas3','sessions');
if ~exist(sessionsRoot, 'dir')
    disp([sessionsRoot ' does not exist. Check that computer is connected to the server']);
    return;
end

folderList = {'D08sessions','D10sessions','D12sessions','D20sessions',...
    'D22sessions','D24sessions','D27sessions'};
numFolders = length(folderList);

% first, combine all the .nex files in each individual subfolder.
% second, move the individual .nex files into a new subfolder to keep them
% out of the way

for iSubject = 7 : numFolders
    % the main directory for each subject
    mainDir = fullfile(sessionsRoot, folderList{iSubject});
    cd(mainDir);
    
    % pull out the subdirectories
    temp = dir;
    numSubdirs = 0;
    if exist('sessionDirs','var')
        clear sessionDirs;
    end
    if length(temp) < 3
        continue;
    end
    for i = 3 : length(temp)
        if strcmp(temp(i).name, '.DS_Store'); continue; end;
        
        if isdir(temp(i).name)
            numSubdirs = numSubdirs + 1;
            sessionDirs{numSubdirs} = temp(i).name;
        end    % if isdir(temp(i))
    end    % for i = 1 : length(temp)
    
    for iSession = 1 : numSubdirs
        curDir = fullfile(mainDir, sessionDirs{iSession});
        finDur = fullfile(curDir, 'Finished');
        cd(finDur);
        nexFiles = dir('*.nex');
        
        if length(nexFiles) > 1
            % if one of the files is a .box.nex file, exclude it
            validNexIdx = [];
            for i = 1 : length(nexFiles)
                if isempty(strfind(nexFiles(i).name, '.box.nex'))
                    validNexIdx = [validNexIdx, i];
                end
            end
            if isempty(validNexIdx)
                disp(['No non-box.nex .nex files in ' curDir, 'skipping...']);
                continue;
            end
            nexFiles = nexFiles(validNexIdx);
            
            if length(nexFiles) > 1
                disp(['more than one .nex file in ' curDir, 'skipping...']);
                continue;
            end
        end
        
        if isempty(nexFiles)
            disp(['no .nex files in ' finDur, 'skipping...']);
            continue;
        end
        
        hsdName = findHSDfile_fromSessionFolder( sessionDirs{iSession} );
        if ~exist(hsdName,'file');continue;end;
        nexName = fullfile(finDur, nexFiles.name);
        if isempty(hsdName)
            disp(['hsd file could not be found for session ' sessionDirs{iSession}]);
            continue;
        end
        DKLextract_wideband(nexName, hsdName);
        
    end    % for iSession = 1 : numSubdirs
        
end    % for iSubject = 1 : numFolders

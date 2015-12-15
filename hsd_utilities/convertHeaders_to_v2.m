function convertHeaders_to_v2( convertDir )

cd(convertDir);

hsdList = dir('*.hsd');
numHSD  = length(hsdList);

for iHSD = 1 : numHSD
    
    hsdName = hsdList(iHSD).name;
    hsdName = fullfile(convertDir, hsdName);
    
    [~, ~] = switchHeader_to_v2( hsdName );
    
end
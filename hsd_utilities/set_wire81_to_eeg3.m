% script to change all wire 81's to eeg3

[fn, pn] = uigetfile({'*.hsd;*.hsdf'}, 'multiselect', 'on');

if ~iscell(fn)
    fn = cellstr(fn);
end

for iFile = 1 : length(fn)
    
    fn{iFile} = fullfile(pn, fn{iFile});
    
    editHSDwires( fn{iFile}, 'name', 81, {'eeg3'});
    
end
function hsdName = findHSDfile_fromDB( c, varargin )
%
% usage: hsdName = findHSDfile_fromDB( c, varargin )
%
% INPUTS:
%   c - single element of a channel structure
%
% varargs:
%
%   'hsddirectory' - directory in which to look for the hsd file
%
% OUTPUTS:
%   hsdName - the full name, including path, of the hsd file. Returns an
%      empty string if the file doesn't exist (at least, can't be found)


% if strcmpi(c.subject(1:3), 'cpc')
%     dirString = ['Greg data' filesep c.subject filesep c.subject '_lfp'];
% else
    switch c.subject
%         case 'IM-149',
%             % data from this rat are ugly, not using them right now
%             dirString = fullfile('/Volumes','dan','Recording data','IM-149_DL-03','IM-149_DL-03_hsd');
        case 'IM-156',
            dirString = fullfile('/Volumes','dan','Recording data','IM-156_DL-10','IM-156_DL-10_hsd');
        case 'IM-157'
            dirString = fullfile('/Volumes','dan','Recording data','IM-157_DL-08','IM-157_DL-08_hsd');
        case 'IM-163',
            dirString = fullfile('/Volumes','dan','Recording data','IM-163_DL-12','IM-163_DL-12_hsd');
        case 'IM-164',
            dirString = fullfile('/Volumes','dan','Recording data','IM-164_DL-22','IM-164_DL-22_hsd');
        case 'IM-166',
            dirString = fullfile('/Volumes','dan','Recording data','IM-166_DL-20','IM-166_DL-20_hsd');
        case 'IM-174',
            dirString = fullfile('/Volumes','dan','Recording data','IM-174_DL-24','IM-174_DL-24_hsd');
        case 'IM-225',
            dirString = fullfile('/Volumes','dan-nas3','Recording data','IM-225_DL-27','IM-225_DL-27_hsd');
        case 'IM-230',
            dirString = ['IM-230_DL-28' filesep 'IM-230_DL-28_lfp'];
    end    % switch c.subject
% end


for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'hsddirectory',
            dirString = varargin{iarg + 1};
    end
end


hsdFile = c.files.highSpeedData.file;
if isempty(hsdFile)
    % file name field was empty
    hsdName = '';
    return;
end
if strcmp(hsdFile(1), '\')
    hsdFile = PCfn2macfn(hsdFile);
end

[~, fn, ext, ~] = fileparts(hsdFile);

fn = fullfile(dirString, [fn ext]);

if exist(fn, 'file')
    hsdName = fn;
else
    hsdName = '';
end
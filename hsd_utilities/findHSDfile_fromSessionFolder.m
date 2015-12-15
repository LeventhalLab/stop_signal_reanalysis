function hsdName = findHSDfile_fromSessionFolder( sessionDir )
%
% usage: hsdName = findHSDfile_fromSessionFolder( sessionDirs )
%
% INPUTS:
%    sessionDir - name of a directory containing session information;
%    should be in the format 'Dxxyyyymmdd', where "D" indicates one of
%    Dan's  rats, xx is the rat number (pre-implant), and yyyymmdd is the
%    date the session was recorded
% varargs:
%
%   'hsddirectory' - directory in which to look for the hsd file
%
% OUTPUTS:
%   hsdName - the full name, including path, of the hsd file. Returns an
%      empty string if the file doesn't exist (at least, can't be found)


switch sessionDir(1:3)
    case 'D10',
        dirString = fullfile('/Volumes','dan','Recording data','IM-156_DL-10','IM-156_DL-10_hsd');
        ratID = 'IM-156';
    case 'D08'
        dirString = fullfile('/Volumes','dan','Recording data','IM-157_DL-08','IM-157_DL-08_hsd');
        ratID = 'IM-157';
    case 'D12',
        dirString = fullfile('/Volumes','dan','Recording data','IM-163_DL-12','IM-163_DL-12_hsd');
        ratID = 'IM-163';
    case 'D22',
        dirString = fullfile('/Volumes','dan','Recording data','IM-164_DL-22','IM-164_DL-22_hsd');
        ratID = 'IM-164';
    case 'D20',
        dirString = fullfile('/Volumes','dan','Recording data','IM-166_DL-20','IM-166_DL-20_hsd');
        ratID = 'IM-166';
    case 'D24',
        dirString = fullfile('/Volumes','dan','Recording data','IM-174_DL-24','IM-174_DL-24_hsd');
        ratID = 'IM-174';
    case 'D27',
        dirString = fullfile('/Volumes','dan-nas3','Recording data','IM-225_DL-27','IM-225_DL-27_hsd');
        ratID = 'IM-225';
end    % switch c.subject


for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'hsddirectory',
            dirString = varargin{iarg + 1};
    end
end

dateNumber = datenum(sessionDir(4:end), 'yyyymmdd');
dateString = datestr(dateNumber, 'yyyy-mm-dd');

testName = fullfile(dirString,[ratID '_' dateString '*.hsd']);

hsdInfo = dir(testName);

if isempty(hsdInfo)
    hsdName = '';
else
    hsdName = fullfile(dirString, hsdInfo.name);
end
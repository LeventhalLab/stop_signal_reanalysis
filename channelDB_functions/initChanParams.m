function cp = initChanParams( varargin )
%
% usage: cp = initChanParams( task, locationName, locationSubclass,
%    subject, date, tetrode, session, channelName )
%
% if no inputs are used, defaults are to task = -1, all other fields
% populated with "any"
%
% default values:
% cp.task = -1;
% cp.locationName = 'any';
% cp.locationSubClass = 'any';
% cp.subject = 'any';
% cp.date = 'any';
% cp.tetrode = 'any';
% cp.session = 'any';
% cp.channelName = 'any';
% cp.isValid = 1;

cp.task = -1;
cp.locationName = 'any';
cp.locationSubClass = 'any';
cp.subject = 'any';
cp.date = 'any';
cp.tetrode = 'any';
cp.session = 'any';
cp.channelName = 'any';
cp.isValid = 1;

cpNames = fieldnames(cp);

for iarg = 1 : nargin
    if ~isempty(varargin{iarg})
        cp.(cpNames{iarg}) = varargin{iarg};
    end
end
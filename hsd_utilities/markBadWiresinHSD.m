function markBadWiresinHSD(fn, varargin)

makePlots = 1;
timeLimits = [100 120];

for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'makeplots',
            makePlots = varargin{iarg + 1};
        case 'timelimits',
            timeLimits = varargin{iarg + 1};
    end
end

hsdHeader = getHSDHeader(fn);

disp('Currently marked bad wires:')
find([hsdHeader.channel.good] == 0)

if makePlots
    plotSessionHSD('filename',fn, 'timelimits', timeLimits);
end

badWires = input('Bad Wires: ');
goodWires = input('Good Wires: ');

editHSDwires(fn, 'good', badWires, zeros(1, length(badWires)));
editHSDwires(fn, 'good', goodWires, ones(1, length(goodWires)));

hsdHeader = getHSDHeader(fn);
disp('good wires:')
[hsdHeader.channel.good]
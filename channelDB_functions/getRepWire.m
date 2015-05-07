function [repWire] = getRepWire( channel )
%
% function to extract the wire number of the representative wire for a
% given channel. If a representative wire has previously been determined,
% that one will be used. If not, the first "good" wire for the channel will
% be used.  Finally, if there are no good wires for this channel, repWire
% is returned as zero.

wires = channel.wire;
repWire = 0;

if isempty(wires.markedGood) && isempty(wires.representativeWire)
    return;
end

if isempty(wires.representativeWire) && max(wires.markedGood) == 0        
    % no good wires for this channel
    return;
end

if isempty(wires.markedGood)
    % no good wires for this channel
    return;
end

if ~isempty(wires.representativeWire)
    if wires.representativeWire==0
    return;
    end
    % representative wire has been marked
    repWire = wires.number(wires.representativeWire);
    return;
end

% make sure a wire marked as bad isn't used (there are some
% of Greg's channel where a wire is marked both "good" and "bad", but
% clearly should be marked "bad")
for iWire = 1 : length(wires.number)
    if isempty(wires.badWires)
        if wires.markedGood(iWire)
            repWire = wires.number(iWire);
            break;
        elseif wires.markedGood(iWire) && ~wires.badWires(iWire)
            repWire = wires.number(iWire);
            break;
        end
    end
end

if strcmpi(channel.tetrode.type,'r')
    if repWire > 2
        % extra check to make sure we don't try to access a wire that isn't
        % available for a reference channel
        repWire = 0;
    end
end

if strcmpi(channel.tetrode.type,'e')
    if repWire > 1
        % extra check to make sure we don't try to access a wire that isn't
        % available for an EEG channel
        repWire = 0;
    end
end


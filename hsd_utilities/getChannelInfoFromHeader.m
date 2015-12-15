function chInfo = getChannelInfoFromHeader( header, wireNum )
%
% function to extract the channel structure information from the hsd header
% for a specific wire
%
% usage: chInfo = getChannelInfoFromHeader( header, wireNum )
%
% INPUTS:
%   header - the hsd header
%   wireNum - the original wire number(s) of interest (likely not the same
%      as the order in the actual file for LFPs)
%
% OUTPUT:
%   chInfo - the channel information for the requested wire. Fields are:
%       .original_number - wire number in the original .hsd file
%       .good - 1 if yes, 0 if no
%       .channel_type - 1 - tetrode, 2 - stereotrode, 3 - ECoG, 4 - probe
%           channel
%       .channel_number - number of this channel among channels of this
%           type (ie, channel_number is 1 for eeg1, even though it's wire
%           27 for a 21-tetrode drive)
%       .bank - which bank the wire is on. For multi-shank probes, use this
%          to say which shank it's on
%       .gain - self-explanatory
%       .low_cut - low frequency in hardware band-pass filter
%       .high_cut - high frequency in hardware band-pass filter
%       .name - the channel's name (ie, eeg1, T01, S1_01, etc)

orig_numbers = [header.channel.original_number];

if ~any(orig_numbers == wireNum)
    chInfo = [];
    disp(['No channels for wire number ' num2str(wireNum)]);
    return
end

chInfo = header.channel(orig_numbers == wireNum);
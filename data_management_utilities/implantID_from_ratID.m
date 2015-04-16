function implantID = implantID_from_ratID(ratID)
%
% usage: 
%
% INPUTS:
%   ratID - a string indicating the rat ID (ie, 'D27', etc)
%
% OUTPUTS:
%   implantID - a string indicating the implant ID corresponding to the
%       rat ID (ie, 'IM225', etc)
ratID_list     = {'D20','D22','D24','D27'};
implantID_list = {'IM166','IM164','IM174','IM225','IM296','IM310','IM328','IM336','IM351','IM365','IM372'};

if strcmpi(ratID(1:2), 'IM')
    implantID = ratID;
else
    ratID_idx = strcmpi(ratID(1:3), ratID_list);

    if isempty(ratID_idx)
        implantID = '';
    else
        implantID = implantID_list{ratID_idx};
    end
end
    
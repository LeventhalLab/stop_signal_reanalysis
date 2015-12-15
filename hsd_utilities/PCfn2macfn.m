function [ mac_filename ] = PCfn2macfn( fn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mac_filename = strrep(fn,'\','/');

if strcmp(mac_filename(1:2), '//')
    mac_filename = strrep(mac_filename,'//','/Volumes/');
end

if ~isempty(strfind(mac_filename, 'neurodata.berkelab.org'))
    mac_filename = strrep(mac_filename,'/neurodata.berkelab.org/','/');
end
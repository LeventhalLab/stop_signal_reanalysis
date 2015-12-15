function writeHSDheader( fn, hsd_header )
%
% usage: writeHSDheader( fn, hsd_header )
%
% function to write an hsd header to a file; useful for replacing a header
% with mistakes in it
%
% INPUTS:
%   fn - filename
%   hsd_header - hsd header structure
%

% % % HSD simple data format for high-speed neurophysiological
% % % data acquisition. Created 2006 by Joshua Berke, Jason Breck, and
% % % Gregory Gage and in the laboratory of Joshua Berke, University
% % % of Michigan, Ann Arbor. Data consists of continuous streams of raw
% % % I16 A/D values, range +-10V, from NI M-series cards via LabVIEW.
% % % 
% % % This text description begins the file.
% % % User-entered description of this experiment starts at 2 * 1024 bytes.
% % % Main header fields begin at 3 * 1024 bytes.
% % % Channel header fields begin at 5 * 1024 bytes.
% % % Data begin at 20 * 1024 bytes.
% % % 
% % % once at the beginning of the HSD file: (begins at 3k)
% % % 
% % %        uint16       Header Version Number (this text describes version 3 !!)
% % %        uint16       Sampling rate (Hz)
% % %        uint16       Number of channels in this file
% % %        char[256] Original name of this file (including the time/date string)
% % %        char[256] Name of auxilliary file 1 (e.g. BOX file)
% % %        char[256] Name of auxilliary file 2 (e.g. AVI file)
% % %        char[256] Name of auxilliary file 3 (e.g. AVI-timestamps file)
% % %        char[256] Name of auxilliary file 4
% % %        char[256] Name of auxilliary file 5
% % %        char[256] Name of auxilliary file 6
% % %        char[100] Subject Name
% % %        char[10]   Date (format "YYYY-MM-DD")
% % %        char[5]    Time (24 hour format "HH:MM")
% % %        char[20]   recordingType - 'tetrode21', 'probe1x32'
% % %        double     downsampled_rate
% % %        char[107] PADDING
% % % 
% % % then, for each channel: (begins at 5k)
% % % 
% % %        uint16      Original Input Wire Number  (e.g. wires 1-81)
% % %        uint16      Good Wire (1 = yes, 0 = No)
% % %        uint16      Channel Type (1=tetrode, 2=stereotrode, 3=eeg, 4=wire)
% % %        uint16      Channel number of this type (0 = none)
% % %        uint16      Wire number of this channel (0 = none)
% % %        uint16      Bank number of this wire (1 = BankA, 2 = BankB, ...)
% % %        uint16      Gain setting for this wire
% % %        uint16      Low cut frequency of hardware filters (Hz)
% % %        uint16      High cut frequency of hardware filters (Hz)
% % %        char[64]  Wire name
% % % 
% % % then, the data (begin at 20k)

% does the file 'fn' already exist?
comment_offset = 2 * 1024;
main_header_offset = 3 * 1024;
channel_header_offset = 5 * 1024;

if ~exist( fn, 'file' )
    % this file does not exist yet; need to write the full header
    fid = fopen(fn, 'w');
    space_pad(1:20480) = ' ';
    fwrite(fid, ['HSD ' space_pad], '*char');
    clear space_pad
    
else
    fid = fopen(fn, 'r+');
    fwrite(fid, 'HSD ', '*char');
end    % end if ~exist( fn )


% write the comment
padLength = main_header_offset - comment_offset;
space_pad(1 : padLength - length(hsd_header.comment)) = ' ';
hComment = [hsd_header.comment space_pad];
fseek(fid, comment_offset, 'bof');
fwrite(fid, hComment, '*char');

% write the main information
fseek(fid, main_header_offset, 'bof');
fwrite(fid, hsd_header.main.version, 'uint16', 0, 'b');
fwrite(fid, hsd_header.main.sampling_rate, 'uint16', 0, 'b');
fwrite(fid, hsd_header.main.num_channels, 'uint16', 0, 'b');

clear space_pad
space_pad(1 : (256 - length(hsd_header.original_filename))) = ' ';
hsd_header.original_filename = [hsd_header.original_filename space_pad];
fwrite(fid, hsd_header.original_filename, '*char', 0, 'b');

for iaux = 1 : 6
    space_pad(1 : (256 - length(hsd_header.aux_filename{iaux}))) = ' ';
    hsd_header.aux_filename{iaux} = [hsd_header.aux_filename{iaux} space_pad];
    fwrite(fid, hsd_header.aux_filename{iaux}, '*char', 0, 'b');
end

clear space_pad
space_pad(1 : (100 - length(hsd_header.subject))) = ' ';
subject_id = [hsd_header.subject space_pad];
fwrite(fid, subject_id, '*char', 0, 'b');

clear space_pad;
space_pad(1 : (10 - length(hsd_header.date))) = ' ';
date_str = [hsd_header.date space_pad];
fwrite(fid, date_str, '*char', 0, 'b');

clear space_pad;
space_pad(1 : (5 - length(hsd_header.date))) = ' ';
time_str = [hsd_header.time space_pad];
fwrite(fid, time_str, '*char', 0, 'b');

if hsd_header.main.version == 3
    clear space_pad;
    space_pad(1 : (20 - length(hsd_header.recordingType))) = ' ';
    recordingType = [hsd_header.recordingType space_pad];
    fwrite(fid, recordingType, '*char', 0, 'b');
end

if isfield(hsd_header.main, 'downsampled_rate')
    fwrite(fid, hsd_header.main.downsampled_rate, 'double', 0, 'b');
else
    fwrite(fid, double(0), 'double', 0, 'b');
end

fseek(fid, channel_header_offset, 'bof');

for iChannel = 1 : hsd_header.main.num_channels

    fwrite(fid, hsd_header.channel(iChannel).original_number, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).good, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).channel_type, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).channel_number, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).wire_number, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).bank, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).gain, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).low_cut, 'uint16', 0, 'b');
    fwrite(fid, hsd_header.channel(iChannel).high_cut, 'uint16', 0, 'b');
    
    clear space_pad;
    space_pad(1 : (64 - length(hsd_header.channel(iChannel).name))) = ' ';
    new_name = [hsd_header.channel(iChannel).name space_pad];
    fwrite(fid, new_name, '*char', 0, 'b');

end


fclose(fid);
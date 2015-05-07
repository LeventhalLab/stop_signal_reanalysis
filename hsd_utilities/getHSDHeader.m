function [ hsd_header ] = getHSDHeader( input_file_name )
% Read New HSD
% 
%  started by Jason, 7/26/06
%  edited by Greg, 14-Nov-2006
%  edited by Dan Leventhal, 9/17/2010

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PARAMETERS                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%samples_per_buffer = 100000;
%   % (this is how many samples to read in and plot)
%channel_to_plot = 21;
%sample_offset = 0; %Fs * time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_header_offset = 3 * 1024;
comment_offset = 2 * 1024;
channel_header_offset = 5 * 1024;
data_offset = 20 * 1024;

if nargin == 0
    [fn pn] = uigetfile('*.hsd', 'Choose a .HSD file');
    input_file_name = [pn fn];
end

dataFilename = input_file_name;

if exist( [input_file_name '.header'], 'file')
    
 
    %Check to see if the data file has a header...
    fin = fopen(input_file_name, 'r');
    magic_string = fread(fin, 4, '*char')';
    data_offset = 0;
    if strcmp(magic_string, 'HSD ')
         data_offset = 20 * 1024;
    end 
    fclose( fin );
    
    input_file_name = [input_file_name '.header'];

end
    
%disp(['  Reading from file ' input_file_name]);

%dataname = [input_file_name ' channel ' num2str(channel) ' starting at ' num2str(sample_offset)];

fin = fopen(input_file_name, 'r');

% Begin reading the file.
%  note that the beginning of the file is the
%  textual description of the file format,
%  and that begins with "HSD "
magic_string = fread(fin, 4, '*char')';
if strcmp(magic_string, 'HSD ')
    %disp('  Magic String "HSD " was found');
else
    disp('  WARNING: Magic String "HSD " not found!');
    disp(['     instead, I found ' magic_string]);
end

% Read file comment 
%  (i.e. the user-entered per-experiment comment field)
fseek(fin, comment_offset, 'bof');
hsd_header.comment = deblank(fread(fin, 1024, '*char')');
%disp('  File comment is as follows:');
%disp(' ');
%disp(hsd_header.comment);
%disp(' ');

% Read in the main HSD file header
fseek(fin, main_header_offset, 'bof');
hsd_header.main.version       = fread(fin, 1, 'uint16', 0, 'b');
hsd_header.main.sampling_rate = fread(fin, 1, 'uint16', 0, 'b');
hsd_header.main.num_channels  = fread(fin, 1, 'uint16', 0, 'b');
%disp(['  I see ' num2str(hsd_header.main.num_channels) ' channels']);
hsd_header.original_filename  = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{1}    = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{2}    = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{3}    = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{4}    = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{5}    = deblank(fread(fin, 256, '*char')');
hsd_header.aux_filename{6}    = deblank(fread(fin, 256, '*char')');
hsd_header.subject  = deblank(fread(fin, 100, '*char')');
hsd_header.date  = deblank(fread(fin, 10, '*char')');
hsd_header.time  = deblank(fread(fin, 5, '*char')');

hsd_header.recordingType = deblank(fread(fin, 20, '*char')');
if hsd_header.main.version ~= 3
    % includes device type
    hsd_header.recordingType = '';
end
hsd_header.main.downsampled_rate = fread(fin, 1, 'double', 0, 'b');

fseek(fin, channel_header_offset, 'bof');
for which_channel = 1:hsd_header.main.num_channels
    
    if hsd_header.main.version == 1
        %        uint16    Original Input Channel Number (before deleting any channels)
        %        uint16    Tetrode number of this channel (0 = none)
        %        uint16    Wire number of this channel (0 = none) 
        %        uint16    Gain setting for this channel
        %        uint16    Low cut frequency of hardware filters (Hz)
        %        uint16    High cut frequency of hardware filters (Hz)
        %        char[64]  Channel Name
        hsd_header.channel(which_channel).original_number = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).good = 1;
        hsd_header.channel(which_channel).channel_type  = 0;
        hsd_header.channel(which_channel).channel_number  = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).wire_number     = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).bank            = 0;
        hsd_header.channel(which_channel).gain            = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).low_cut         = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).high_cut        = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).name            = deblank(fread(fin, 64, '*char')');
       
    else %assume version 2
        hsd_header.channel(which_channel).original_number = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).good  = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).channel_type  = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).channel_number  = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).wire_number     = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).bank            = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).gain            = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).low_cut         = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).high_cut        = fread(fin, 1, 'uint16', 0, 'b');
        hsd_header.channel(which_channel).name            = deblank(fread(fin, 64, '*char')');
    end
end

hsd_header.dataOffset = data_offset;
hsd_header.dataFilename = dataFilename;

% Now, read in a chunk of data, and plot it just for... err... fun...
% fseek(fin, data_offset + 2 * sample_offset * hsd_header.main.num_channels, 'bof');
%[input_buffer, number_elements_read] = fread(fin, [hsd_header.main.num_channels, samples_per_buffer], 'int16', 0, 'b');

fclose(fin);

%plot(input_buffer(channel_to_plot));
%xlim([0 200]);


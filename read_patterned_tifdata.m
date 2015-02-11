function [y, s] = read_patterned_tifdata(s, frames)
%read_patterned_tifdata
%Read data from a tif file with a repeating structure, by guessing the location of each frame's data.
%This function only works with multipage tifs that have the same amount of data for each image,
%stored in the same way and with the same between-frames-offset for each each pair of consecutive frames.
%Because this function reads the necessary data for each frame directly from the file without reading
%offset or tags for each frame, it can be much faster than standard tif readers.
%
%simplest use:
%data = read_patterned_tifdata(filename)
%reads all data from the file into an array index as
%data(row, column, frame, channel)
%
%data = read_patterned_tifdata(filename, frames)
%specifies which frames to load
%
%[data,s] = read_patterned_tifdata(filename, frames)
%also returns a structure containing an open file handle (for used with fread etc.)
%Afterwards one can use:
%data = read_patterned_tifdata(s, frames) to read data without parsing/scanning the file,
%this can give a significant speedup.
% ---> When finished, the user should call fclose(s.fid) !!! <---
%If a call to read is made with a filename instead of a structure as the first argument,
%then the file handle is left open only if a second output is returned.
%If a call that should return a file handle fails (due to bad filename, invalid file, etc.)
%then there should not be a file handle left open.
%
%s also contains many other fields providing information about the tif file and where to find its image data.
%For large files, it's best to initialize without reading any data:
%[~,s] = read_patterned_tifdata(filename, [])
%
%Planar and interleaved tif configurations are both supported, but compression is not.
%
%The first 3 frames of the file are anlayzed in detail to verify that a repeating pattern is present.
%However, there is no absolute guarantee that the returned data will be corrected if the file is
%patterned on the first 3 frames and then departs from that pattern afterwards
%
%This function works on abberant tif files larger than 4GB (standard tif format, not bigtiff files!).
%For example, as of 29.01.2015, the beta version of the Scanimage 5 microscope software (http://vidriotechnologies.com/)
%will save such files if instructed to record more data than can fit within the standard 4GB tif limit.
%Obviously, such files violate the tiff standard and must contain invalid values for some 32-bit fields
%that store offsets into the file itself. However, since this function relies on predicting data location
%instead of reading it, it can read data from these files without issues.
%
%This function makes use of ideas and some code from tiffread, version 2.4:
%Francois Nedelec, nedelec (at) embl.de,
%Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
%http://www.embl.org
%http://www.cytosim.org
%
%Written by David Greenberg January 2015
%david.greenberg@caesar.de
filehandlecreated = false;
if isa(s, 'char')
    s = open_patterned_tifdata(s);
    filehandlecreated = true;
elseif ~s.file_is_open && exist(s.filename, 'file') == 2 %this lets us save s, restart the computer, load s, and call a read function without dealing with tags etc.
    s.fid = fopen(s.filename, 'r'); %this assumes that the user hasn't modified the tif file. but if the user isn't sure they can always pass the filename instead of the struct.
    filehandlecreated = true;
end
if nargin < 2
    frames = 1:s.nframes; %may cause out of memory errors on large files
else
    assert(all(frames <= s.nframes), 'frame index exceeds number of images in tif file');
end
if filehandlecreated
    try
        y = read_patterned_tifdata_readframes(s, frames); %no try/throw/catch, hopefully this gives a speed boost
    catch ex
        fclose(s.fid); %clean up our mess
        ex.rethrow();
    end
else
    y = read_patterned_tifdata_readframes(s, frames); %no try/throw/catch, hopefully this gives a speed boost
end
if filehandlecreated && nargout < 2
    fclose(s.fid);
    s.file_is_open = false;
end

function y = read_patterned_tifdata_readframes(s, frames)
y = zeros([s.imagesize(1:2) numel(frames) s.imagesize(3)], s.data_matlab_class);
nstrips = numel(s.tifinfo.StripByteCounts);
npixels = prod(s.imagesize(1:2));
csbytes = cumsum([0 s.tifinfo.StripByteCounts]);
cssamples = csbytes / s.bytespersample;
df = diff(frames);

if nstrips == 1 && numel(frames) > 1 && df(1) > 0 && all(df == df(1))
    %in this case we are making evenly sized file reads of the same size, so we can consolidate all file freads into a single fread command
    samplesperframe = prod(s.imagesize);
    bytes2skip = df(1) * s.offset_bytes_per_frame - samplesperframe * s.bytespersample; %bytes between the end of one frame and the beginning of the next    
    so = s.tifinfo.StripOffsets + s.offset_bytes_per_frame * (frames(1) - 1); %offset for the first and only strip of the first image to be read
    assert(fseek(s.fid, double(so), -1) == 0, 'seek failed');
    nsamples2read = samplesperframe * numel(frames);
    precision = sprintf('%d*%s=>%s', samplesperframe, s.data_matlab_class, s.data_matlab_class);
    dataraw_all = fread(s.fid, nsamples2read, precision, double(bytes2skip), s.machineformat); %read all frames and channels in a single go
    %now we have to reorder the data to use matlab indexing etc.
    for u = 1:numel(frames) %FIXME it might not be crazy to do this loop with mex/openMP
        frameoffset = (u - 1) * samplesperframe;
        for c = 1:s.imagesize(3)            
            if s.imagesize(3) == 1 || s.tifinfo.PlanarConfiguration == 2
                channeloffset = frameoffset + (c - 1) * npixels;
                y(:, :, u, c) = reshape( ...
                    dataraw_all(channeloffset + 1:channeloffset + npixels), ...
                    [s.imagesize(2) s.imagesize(1)])'; %we need the transpose since tifs store by row first and matlab stores by column first                
            else
                y(:, :, u, c) = reshape( ...
                    dataraw_all(frameoffset + (0:npixels - 1) * s.imagesize(3) + c), ...
                    [s.imagesize(2) s.imagesize(1)])'; %we need the transpose since tifs store by row first and matlab stores by column first                
            end            
        end
    end
else %general cases, requiring one call to fread for each strip of each frame
    dataraw = zeros([prod(s.imagesize) 1], s.data_matlab_class); %will store raw data array for one image, of the correct numeric type but not initially in the correct order
    for u = 1:numel(frames)
        so = s.tifinfo.StripOffsets + s.offset_bytes_per_frame * (frames(u) - 1); %offsets for each strip of image ii(u)
        for jj = 1:nstrips
            assert(fseek(s.fid, double(so(jj)), -1) == 0, 'seek failed');
            nsamples2read = s.tifinfo.StripByteCounts(jj) / s.bytespersample;
            [dataraw(cssamples(jj) + 1:cssamples(jj + 1)), nsamplesread] = fread(s.fid, nsamples2read, s.data_matlab_class, s.machineformat);
            assert(nsamplesread == nsamples2read, 'failed to read all expected data from tif file');
        end
        for c = 1:s.imagesize(3)
            if s.imagesize(3) == 1 || s.tifinfo.PlanarConfiguration == 2
                dataraw_nextchannel = dataraw((c - 1) * npixels + 1:c * npixels);
            else
                dataraw_nextchannel = dataraw(c:s.imagesize(3):end);
            end
            y(:, :, u, c) = reshape(dataraw_nextchannel, [s.imagesize(2) s.imagesize(1)])'; %we need the transpose since tifs store by row first and matlab stores by column first
        end
    end
end

function offset_bytes_per_frame = verify_repeating_data_pattern(tifinfo, tifobj)
so2 = tifobj.getTag('StripOffsets');
bc2 = tifobj.getTag('StripByteCounts');
assert(numel(so2) == numel(tifinfo.StripOffsets) && numel(bc2) == numel(tifinfo.StripByteCounts), 'images must have the same number of strips');
assert(all(so2 > tifinfo.StripOffsets), 'images cannot overlap');
assert(all(bc2 == tifinfo.StripByteCounts), 'images must have the same byte count per strip');
d = so2 - tifinfo.StripOffsets;
assert(all(d == d(1)), 'image strip offsets must have the same spacing across images');
offset_bytes_per_frame = uint64(d(1));
%find the maximum number of frames that fit in the file

function s = open_patterned_tifdata(filename)
assert(exist(filename,'file') == 2, 'file not found');
fileinfo = dir(filename);
filesize = uint64(fileinfo.bytes);
tifobj = Tiff(filename);
tifinfo = get_tifinfo(tifobj);
assert(all(tifinfo.StripOffsets + tifinfo.StripByteCounts) < filesize, 'file is too small for the offsets it contains');
if tifobj.lastDirectory() %file contains only one frame
    offset_bytes_per_frame = 0;
    nframes = 1;
else %advance to the 2nd directory and check whether the data spacing is patterned
    tifobj.nextDirectory(); 
    offset_bytes_per_frame = verify_repeating_data_pattern(tifinfo, tifobj);
    if ~tifobj.lastDirectory() %try a 3rd directory too if the file contains at least 3 directories
        tifobj.nextDirectory();
        offset_bytes_per_frame2 = verify_repeating_data_pattern(tifinfo, tifobj); %this offset is the distance between frames 1 and 3
        assert(offset_bytes_per_frame * 2 == offset_bytes_per_frame2, 'offset between frames 1-2 is not the same as offset between frames 2-3');
    end    
    bytesreadfromoffset = max(tifinfo.StripOffsets + tifinfo.StripByteCounts); %to read frame n, the file size must be at least (n * offset_bytes_per_frame + bytesreadfromoffset) bytes
    nframes = double(1 + floor(double(filesize - bytesreadfromoffset) / double(offset_bytes_per_frame)));    
end
[tifinfo.StripOffsets, tifinfo.StripByteCounts] = consolidate_strips(tifinfo.StripOffsets, tifinfo.StripByteCounts);
fid = fopen(filename, 'r'); %FIXME memmapfile might give better performance?
assert(fid > 0, 'failed to open file');
try
    s = orderfields(struct( ....
        'filename',                 which(filename), ... %try to get the full path
        'filesize',                 filesize, ...
        'tifinfo',                  tifinfo, ...
        'fid',                      fid, ...
        'imagesize',                double([tifinfo.ImageLength tifinfo.ImageWidth tifinfo.SamplesPerPixel]), ...
        'machineformat',            patterned_tif_getmachineformat(fid), ...
        'file_is_open',             true, ...
        'offset_bytes_per_frame',   offset_bytes_per_frame, ...
        'nframes',                  nframes, ...
        'bytespersample',           tifinfo.BitsPerSample / 8, ...
        'data_matlab_class',        patterned_tif_getdataclass(tifinfo) ...
        ));
catch ex
    fclose(fid); %clean up our mess
    ex.rethrow();
end

function tifinfo = get_tifinfo(tifobj)
tn = tifobj.getTagNames;
tifinfo = struct();
for u = 1:numel(tn)
    try %#ok<TRYNC>
        v = tifobj.getTag(tn{u});
        tifinfo.(tn{u}) = v;
    end
end
tifinfo.StripOffsets = reshape(tifinfo.StripOffsets, 1, []); %just in case matlab's Tiff implementation changes this in the future
tifinfo.StripByteCounts = reshape(tifinfo.StripByteCounts, 1, []); %just in case matlab's Tiff implementation changes this in the future
tifinfo = orderfields(tifinfo);
check_patterned_tifinfo(tifinfo);

function [so, bc] = consolidate_strips(so, bc)
%combine directly adjacent strips to reduce the number of calls to fread when retrieving data
%idea for doing this step comes from tiffread2 by Francois Nedelec, though the implementation below which avoids loops is new
mergewithprev = [false, so(1:end - 1) + bc(1:end - 1) == so(2:end)];
totalbytes = sum(bc);
bc_exclusiveprefixsum = cumsum([0 bc(1:end - 1)]); %total number of bytes BEFORE each strip, not counting the strip itself
%merge strips:
so(mergewithprev) = [];
bc_exclusiveprefixsum(mergewithprev) = [];
bc = diff([bc_exclusiveprefixsum totalbytes]);

function dataclass = patterned_tif_getdataclass(tifinfo)
switch tifinfo.SampleFormat
    case 1
        dataclass = sprintf('uint%i', tifinfo.BitsPerSample);
    case 2
        dataclass = sprintf('int%i', tifinfo.BitsPerSample);
    case 3
        if tifinfo.BitsPerSample == 32
            dataclass = 'single';
        else
            assert(tifinfo.BitsPerSample == 64, 'BitsPerSample must be 32 or 64 for SampleFormat 3 (floating point)');
            dataclass = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', TIF.SampleFormat);
end

function machineformat = patterned_tif_getmachineformat(fid)
frewind(fid);
%this function is based mostly on tiffread2 by Francois Nedelec
[byte_order, nread] = fread(fid, 2, '*char');
assert(nread == 2, 'failed to read machine format');
if ( strcmp(byte_order', 'II') )
    machineformat = 'l';                                %normal PC format
elseif ( strcmp(byte_order','MM') )
    machineformat = 'b';
else
    error('This is not a TIFF file (no MM or II).');
end

function check_patterned_tifinfo(tifinfo) %valid the tif file as making sense and being a file we can use
mandatory_fields = {'StripOffsets', 'StripByteCounts', 'BitsPerSample', 'ImageWidth', 'ImageLength', 'SampleFormat', 'Compression', 'SamplesPerPixel'};
assert(all(ismember(mandatory_fields, fieldnames(tifinfo))), 'one or more required tags were not found in the first directory of the tif file');
assert(tifinfo.Compression == 1, 'compressed data is not supported');
if tifinfo.Compression == 1
    assert(tifinfo.ImageWidth * tifinfo.ImageLength * tifinfo.SamplesPerPixel * (tifinfo.BitsPerSample / 8) == sum(tifinfo.StripByteCounts), 'incorrect amount of data per image');
end
if tifinfo.SamplesPerPixel ~= 1
    assert(isfield(tifinfo,'PlanarConfiguration'), 'Planar configuration must be specified when SamplesPerPixel > 1'); %Nedelec in tiffread2: "Contributed by Stephen Lang"
end
assert(numel(tifinfo.StripOffsets) > 0, 'no image data in this tif file');
assert(numel(tifinfo.StripOffsets) == numel(tifinfo.StripByteCounts),'invalid configuration');
if numel(tifinfo.BitsPerSample) > 1
    assert(all(tifinfo.BitsPerSample == tifinfo.BitsPerSample(1)), 'BitsPerSample must be the same for all image channels');
    tifinfo.BitsPerSample = tifinfo.BitsPerSample(1);
end
assert(all(diff(tifinfo.StripOffsets) > 0), 'StripOffsets must be an increasing sequence'); %note that all([]) == true so this works when there's only one offset
assert(double(tifinfo.BitsPerSample) / 8 == round(double(tifinfo.BitsPerSample) / 8),'bits per sample must be a multiple of 8');
samplesperstrip = double(tifinfo.StripByteCounts) / (double(tifinfo.BitsPerSample) / 8);
assert(all(samplesperstrip == round(samplesperstrip)), 'the size in bytes of each strip must be an integer number of the number of bytes per sample');

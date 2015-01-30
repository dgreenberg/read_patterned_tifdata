# read_patterned_tifdata
A matlab function for rapidly reading multipage tif files with a repeating structure by guessing where the data will be. Can read standard-violating read files larger than 4GB.

#Help comments from the .m file:
read_patterned_tifdata

Read data from a tif file with a repeating structure, by guessing the location of each frame's data.
This function only works with multipage tifs that have the same amount of data for each image,
stored in the same way and with the same between-frames-offset for each each pair of consecutive frames.
Because this function reads the necessary data for each frame directly from the file without reading
offset or tags for each frame, it can be much faster than standard tif readers.

simplest use:

data = read_patterned_tifdata(filename)

reads all data from the file into an array index as
data(row, column, frame, channel)

data = read_patterned_tifdata(filename, frames)

specifies which frames to load

[data,s] = read_patterned_tifdata(filename, frames)

also returns a structure containing an open file handle (for used with fread etc.)
Afterwards one can use:

data = read_patterned_tifdata(s, frames) to read data without parsing/scanning the file,

this can give a significant speedup.

---> When finished, the user should call fclose(s.fid) !!! <---
If a call to read is made with a filename instead of a structure as the first argument,
then the file handle is left open only if a second output is returned.
If a call that should return a file handle fails (due to bad filename, invalid file, etc.)
then there should not be a file handle left open.

s also contains many other fields providing information about the tif file and where to find its image data.
For large files, it's best to initialize without reading any data:

[~,s] = read_patterned_tifdata(filename, [])

Planar and interleaved tif configurations are both supported, but compression is not.

The first 3 frames of the file are anlayzed in detail to verify that a repeating pattern is present.
However, there is no absolute guarantee that the returned data will be corrected if the file is
patterned on the first 3 frames and then departs from that pattern afterwards

This function works on abberant tif files larger than 4GB (standard tif format, not bigtiff files!).
For example, as of 29.01.2015, the beta version of the Scanimage 5 microscope software (http://vidriotechnologies.com/)
will save such files if instructed to record more data than can fit within the standard 4GB tif limit.
Obviously, such files violate the tiff standard and must contain invalid values for some 32-bit fields
that store offsets into the file itself. However, since this function relies on predicting data location
instead of reading it, it can read data from these files without issues.

This function makes use of ideas and some code from tiffread, version 2.4:
Francois Nedelec, nedelec (at) embl.de,
Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
http://www.embl.org
http://www.cytosim.org

Written by David Greenberg January 2015
david.greenberg@caesar.de

function [oldHeader, newHeader] = switchHeader_to_v2( fn )

oldHeader = getHSDHeader( fn );

newHeader = oldHeader;
newHeader.main.version = 2;

writeHSDheader_ver2(fn, newHeader);

newHeader = getHSDHeader(fn);
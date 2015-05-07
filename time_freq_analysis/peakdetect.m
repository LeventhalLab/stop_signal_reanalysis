function [peakLocations, peakValues] = peakDetect(x)

% function to detect peaks in a signal

dx = diff(x) ./ abs(diff(x));

ddx = diff(dx);

peakLocations = find(ddx == -2) + 1;

peakValues = x(peakLocations);
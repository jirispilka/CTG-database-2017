function [data, data_begin_end, seg_remainder] = extractSegments(x,winlen_samp,overlap_samp)
% Extract segmets of predifined length and overlap
%
% Inputs:
%  x             [Nx1] - signal
%  winlen_samp   [int] - window length in samples
%  overlap_samp  [int] - step in samples
%
% Outputs:
%  data             [LxD] - output segments,
%                         - L length of segments, D - number of segments
%  data_begin_end   [Dx2] - markes of begins and ends of segments
%                           - D - number of segments, 2 columns [begin, end]
%
% Jiri Spilka, 11/2013
% ENS Lyon

% if a windows is bigger than data length than use all data
if isinf(winlen_samp) || winlen_samp > length(x)
    winlen_samp = length(x);
    overlap_samp = 0;
end

% number of segments
k = (length(x) - overlap_samp)/(winlen_samp - overlap_samp);

% Uncomment the following line to produce a warning each time the data
% segmentation does not produce an integer number of segments.
% if fix(k) ~= k
%   warning('extractSegments:MustBeInteger','The number of segments is not an integer, truncating data.');
% end

fRemainder = k - fix(k);
k = fix(k);

% init arrays
data = zeros(winlen_samp,k);
LminusOverlap = winlen_samp - overlap_samp;
xstart = 1:LminusOverlap:k*LminusOverlap;
xend = xstart+winlen_samp-1;

seg_remainder = [];
data_begin_end = [xstart', xend'];

for i = 1:k
    data(:,i) = x(xstart(i):xend(i));
    
    if i == k && fRemainder ~= 0
        seg_remainder = x(xend(i):end);
    end
end
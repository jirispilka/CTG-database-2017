function [aFhrInterp] = interpolateSmallGaps(aFhr,nFs, bPlot, nSmallGapLim_samp)
% INTERPOLATESMALLGAPS interpolate gaps shorter then 15 seconds
%
% Synopsis:
%  [record] = interpolateSmallGaps(aFhr [,nFs [,bPlot]])
%
% Description:
%  Interpolate gaps shorter then 15 seconds with spline interpolation
%  (pchip).
%  
% Input:
%  aFhr - [nx1] FHR signal
%  nFs - [int][optional] sampling frequency [Hz]
%  bPlot - [bool][optional] debug 
%
% Output:
%  aFhrInterp - [nxdouble] interpolated FHR signal
%
% Example:
%
% See also:
%
% About:
% Jiri Spilka
% Czech Technical University in Prague, 2015
%
% Modifications:
%

if ~exist('nFs','var')
    nFs = 4;
end

if ~exist('bPlot','var')
    bPlot = 0;
end

if ~exist('nSmallGapLim_samp','var')
    nSmallGapLim_samp = nFs*5;
end

%% find small gaps

[aGapStart aGapLength] = findFhrGaps(aFhr);

idxSmall = find(aGapLength <= nSmallGapLim_samp);
idxNotSmall = find(aGapLength > nSmallGapLim_samp);

% if a signal begins or ends with NaNs --> do not interpolate there
[~,aGapAtBegin, aGapAtEnd] = removeNaNsAtBeginAndEnd(aFhr);

if bPlot
    figure(10350); clf
    plot(aFhr)
    hold on;
    stem(aGapStart(idxSmall), aGapLength(idxSmall),'k');
    hold on;
    stem(aGapStart(idxNotSmall), aGapLength(idxNotSmall),'r');
    grid on;
    legend('FHR','small gap','large gap')
end

% interpolate
idxSamples = ~isnan(aFhr);
aFhrWithouNans = aFhr(idxSamples);
xTimeOld = find(idxSamples);
xTimeNew = 1:length(aFhr);
aFhrInterp = pchip(xTimeOld,aFhrWithouNans,xTimeNew);

% clean begin and end of signal
aFhrInterp(aGapAtBegin) = NaN;
aFhrInterp(aGapAtEnd) = NaN;

% delete big gaps
for i = idxNotSmall'
    nFrom = aGapStart(i);
    nTo = aGapStart(i) + aGapLength(i);
    aFhrInterp(nFrom:nTo) = NaN;
end

if size(aFhrInterp,2) ~= 1
    aFhrInterp = aFhrInterp';
end

% if bPlot
%     figure
%     plot(aFhrInterp,'r');
%     hold on
%     plot(aFhr);
% end

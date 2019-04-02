function [aFhrInterp,aGapAtBegin,aGapAtEnd] = interpolateAllGaps(aFhr, nFs, bPlot)
% INTERPOLATEALLGAPS interpolate all missing data
%
% Synopsis:
%  aFhrInterp = interpolateAllGaps(aFhr, bPlot)
%
% Description:
%  Interpolate all missing data using piecewise cubic 
%  interpolation. (pchip - Matlab builtin function)
%
% Input:
%  aFhr   - [nx1 int] input data
%  nFs    - [int] sampling freqeuncy
%  bPlot  - debug plotting
%
% Output:
%  aFhrInterp  - [nx1 int] interpolated data
%
% Examples:
%
% See also:
%  INTERPOLATESMALLGAPS
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

% if a signal begins or ends with NaNs --> do not interpolate there
[~, aGapAtBegin, aGapAtEnd] = removeNaNsAtBeginAndEnd(aFhr);

% interpolate
idxNotNaN = ~isnan(aFhr);
aFhrWithouNans = aFhr(idxNotNaN);

xTimeOld = find(idxNotNaN);
xTimeNew = (1:length(aFhr))';
%xTimeNew = (1:1/nFs:length(aFhr))';
aFhrInterp = pchip(xTimeOld,aFhrWithouNans,xTimeNew);

if nFs == 4
    % clean begin and end of signal
    aFhrInterp(aGapAtBegin) = NaN;
    aFhrInterp(aGapAtEnd) = NaN;
else
    if ~isempty(aGapAtBegin)
        d = (xTimeNew >= aGapAtBegin(1) & xTimeNew <= aGapAtBegin(end));
        aFhrInterp(d) = NaN;
    end
    
    if ~isempty(aGapAtEnd)
        d = (xTimeNew >= aGapAtEnd(1) & xTimeNew <= aGapAtEnd(end));
        aFhrInterp(d) = NaN;
    end    
end

if bPlot
    figure
    plot(xTimeNew,aFhrInterp,'xr--');
    hold on
    plot(xTimeOld,aFhrWithouNans,'.b--');
end

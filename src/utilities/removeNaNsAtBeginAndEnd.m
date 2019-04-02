function [aFhrOut, aGapAtBegin, aGapAtEnd] = removeNaNsAtBeginAndEnd(aFhr)
% remove NaNs at bebin and end of signal
%
% Jiri Spilka
% ENS Lyon, 2014

aFhrOut = aFhr;

% if a signal begins or ends with NaNs --> do not interpolate there
aGapAtBegin = [];
idxSamples = find(~isnan(aFhr));
if ~isempty(idxSamples)
    if idxSamples(1) ~= 1 % if signal not start at 
        aGapAtBegin = 1:idxSamples(1)-1;
    end
end

idxSamples = find(~isnan(aFhr));
aGapAtEnd = [];
if ~isempty(idxSamples)
    if idxSamples(end) ~= length(aFhr)
        aGapAtEnd = idxSamples(end)+1:length(aFhr);
    end
end

% remove gap at the begin
% remove gap at the end
aFhrOut([aGapAtBegin, aGapAtEnd]) = [];
function [aGapStart, aGapLength] = findFhrGaps(aFhrIn)
% FINDFHRGAPS find all gaps in a FHR signal
%
% Synopsis:
%  [aGapStart aGapLength] = findFhrGaps(aFhrIn)
%
% Description:
%  Find all gaps in a FHR signal
%
% Input:
%  aFhrIn - [nx1] FHR signal
%
% Output:
%  aGapStart - [nx1] begin of gaps
%  aGapLength - [nx1] 
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
aGapStart = zeros(length(aFhrIn),1);
aGapLength = zeros(length(aFhrIn),1);
tempGapStart = 0;
tempGapLength = 0;

ncnt = 0;
ncnt_nogap = 0;

for i = 1:length(aFhrIn)
    % if there is a gap
    if aFhrIn(i) == 0 || aFhrIn(i) == -1 || isnan(aFhrIn(i)) || isinf(aFhrIn(i))
        if tempGapStart == 0
            tempGapStart = i; % save the first index
            ncnt = ncnt + 1;
        end
        tempGapLength = tempGapLength + 1;
    else
        % if there was a gap and now the sample is normal
        if tempGapStart ~= 0
            aGapStart(ncnt) = tempGapStart;
            aGapLength(ncnt) = tempGapLength;
            tempGapStart = 0;
            tempGapLength = 0;
        else
            ncnt_nogap = ncnt_nogap + 1;
        end
    end
end

% the gap at the end of signal is not really gap -> discard
if tempGapStart ~= 0
    ncnt = ncnt - 1;
end

aGapStart = aGapStart(1:ncnt,1);
aGapLength = aGapLength(1:ncnt,1);
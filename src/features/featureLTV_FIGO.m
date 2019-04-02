function [LTV_Delta_total] = featureLTV_FIGO(inDataBpm,nFs)
% extract long term variability/irregularity (LTV) based on FIGO
%
% Synopsis:
%  [LTV_Delta_total] = featureLTV_FIGO(inDataBpm,nFs)
%
% Description:
%
% Input:
%  inDataBpm - [nx1] vector of FHR values in bpm
%  nFs - [int] sampling frequencyvariability 
%
% Output:
%  cFeatures - [struct] values of extracted features
%
% Examples:
%
% See also: 
%
% About: 
% Jiri Spilka
% Czech Technical University in Prague, 2014
%
% Modifications:
%  j.s.  10.1.2013 - renamed to delat total
%  v.ch. 01.10.2012 - added output arguments and output in bpm

nLength_60s = 60*nFs;

% Divide input FHR to 60 seconds segments. 
nNrSubIntervals = fix(length(inDataBpm)/nLength_60s);
aTemp_bpm = zeros(nNrSubIntervals,1);

cnt = 0;
for i = 1:nNrSubIntervals
    nBegin = (i-1)*nLength_60s+1;
    nEnd = (i)*nLength_60s;
    inDataTempBpm = inDataBpm(nBegin:nEnd);
    
    if missingMoreThan50percent(inDataTempBpm)
        continue;
    end
    
    cnt = cnt + 1;
    %% Dawes-LTV Long-term variability
    aTemp_bpm(cnt) = nanmax(inDataTempBpm) - nanmin(inDataTempBpm);
end

% remove alocated space that was not used
aTemp_bpm = aTemp_bpm(1:cnt);

LTV_Delta_total = round(nanmean(aTemp_bpm));
%LTV_Dawes_bpm_min = round(min(aDawesTemp_bpm));
%LTV_Dawes_bpm_max = round(max(aDawesTemp_bpm));

%%
function b = missingMoreThan50percent(inData)

b = false;
nLengthMissing = length(find(isnan(inData)));
if nLengthMissing > length(inData)/2
    b = true;
end
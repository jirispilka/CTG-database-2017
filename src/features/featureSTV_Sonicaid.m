function [STV_Sonicaid_ms] = featureSTV_Sonicaid(inDataBpm, nFs, sFhrType)
% FEATURESTV_Sonicaid extracts short term variability (STV) features
% according to Sonicaid formula: [Reference needed here]
%
% Synopsis:
%  STV_Sonicaid_ms = featureSTV_Sonicaid(inDataBpm , nFs , sFhrType])
%
% Description:
%  Algorithm implemented according to Pardey 2002
%   J. Pardey, M. Moulden, and C. W. G. Redman. A computer system for the numerical analysis of nonstress tests.
%   Am J Obstet Gynecol, 186(5):1095–1103, May 2002.
%
% Input:
%  inDataBpm - [nx1] vector of FHR values in bpm
%  nFs - [int] sampling frequency
%  sFhrType - [string][optional] type of signal: 'US','FECG'
%
% Output:
%  STV_Sonicaid_ms - [float] value in ms
%
% Examples:
%
%
% See also:
%  TESTFEATURESSTV
%
% About:
% Jiri Spilka
% Czech Technical University in Prague, 2015
%
% Modifications:
%

% if ~exist('sFhrType','var')
%      sFhrType = 'US';
% end

inDataRR_original = convertBpmRR(inDataBpm);
nNrSamples = length(inDataRR_original);
nLength_60s_samp = nFs*60;

if nNrSamples < nLength_60s_samp
    error('Not enough data for STV computation, minimal length is 60sec');
end

%% Sonicaid
nNrSubIntervals = fix(length(inDataRR_original)/nLength_60s_samp);
aSonicAidTemp = zeros(nNrSubIntervals,1);

bEvenlySampled = true;

cnt = 0;
for i = 1:nNrSubIntervals
    nBegin = (i-1)*nLength_60s_samp+1;
    nEnd = (i)*nLength_60s_samp;
    dataTempRR = inDataRR_original(nBegin:nEnd);
    
    if missingMoreThan50percent(dataTempRR)
        continue;
    end
    
    %if strcmp(sFhrType,'US')
        tempStv = sonicaidStv(dataTempRR, nFs, bEvenlySampled);
    %else
    %    tempStv = sum(abs(diff(dataTempRR)))/length(dataTempRR);
    %end
    
    if ~isnan(tempStv)
        cnt = cnt + 1;
        aSonicAidTemp(cnt) = tempStv;
    end
end

STV_Sonicaid_ms = mean(aSonicAidTemp(1:cnt));
%acMinuteFeat.Sonicaid = aSonicAidTemp(1:cnt);

%% function to compute sonicaid stv
function fStv = sonicaidStv(inDataRR, nFs, bEvenlySampled)
% compute stv based on sonicaid

global fAverageOverSampSec;
fAverageOverSampSec = 3.75;

% init
h = 1; % count subintervals
nNrSamples = length(inDataRR);
nNrSubIntervals = nNrSamples;

Rm = zeros(nNrSubIntervals,1);

nBegin = h;
% r = number of RR intervals in 3.75 s - value of 3.75 because of easy
% computation - 60/16 = 3.75sec.
nEnd = findEnd(inDataRR, nBegin, nFs, bEvenlySampled);

while nEnd <= nNrSamples
    
    aTempData = inDataRR(nBegin:nEnd);
    % there is a mistake in the referenced article. For evenly spaced FHR
    % series, the m should be equal to m = mean(FHR) in 3.75s but this
    % number is far more larger than m = number of beats for unevenly
    % spaced data.
    %m = length(aTempData);
    
    %     Rm(h) = sum(abs(diff(aTempData))./m);
    Rm(h) = mean(aTempData);
    
    nBegin = nEnd + 1;
    
    % if the end of data was reached -> terminate
    if nEnd < nNrSamples
        h = h + 1;
        nEnd = findEnd(inDataRR, nBegin,nFs, bEvenlySampled);
    else
        %aTempData = inDataRR(nBegin:nEnd);
        %Rm(h) = sum(abs(diff(aTempData))./m);
        nEnd = nEnd + 1;
    end
end
Rm = Rm(1:h);
Rm = Rm(~isnan(Rm)); % remove those that contained NaN values
Rm = Rm(~(Rm == 0)); % remove those that were not computed 
% (typically not enough samples at the end of a record)
h = length(Rm);

if h ~= 0 
    fStv = sum(abs(diff(Rm)))/h;
else % if are all NaN
    fStv = NaN;
end

%% find range of computation
function nEnd = findEnd(inDataRR,h,nFs,bEvenlySampled)

global fAverageOverSampSec;

if bEvenlySampled
    r = fAverageOverSampSec*nFs;
    nEnd = h + round(r) - 1;
else
    idxEnd = h - 1 + find(cumsum(inDataRR(h:end)) < fAverageOverSampSec*1000);
    if isempty(idxEnd)
        idxEnd = length(inDataRR) + 1;
    end
    nEnd = idxEnd(end);
end

%%
function b = missingMoreThan50percent(inData)

b = false;
nLengthMissing = length(find(isnan(inData)));
if nLengthMissing > length(inData)/2
    b = true;
end
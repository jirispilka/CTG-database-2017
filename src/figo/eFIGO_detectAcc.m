function [accVector,accPos] = eFIGO_detectAcc(nxData, baseline, nFs, ampTreshold,lengthTreshold, bPlot)

% [diagBaaseline_out, decelBaseline_out] = eFIGO_detectAcc(fhr,fs,bPlotAlgoSteps)
%
% Detection of FHR baseline
%
% INPUTS:
%       fhr - vector of size (X,1)
%       fs - sampling frequency in Hz
%       bPlotAlgoSteps - fully illustrated passage through the code
%
% OUTPUTS:
%       diagBaaseline_out - baseline usable for diagnostics
%       decelBaseline_out - baseline usable for detection of decelerations/acceletrations%
%
% Authors:
%       Jiri Spilka, Vaclav Chudacek
%       http://ctg.ciirc.cvut.cz
%       CIIRC, CVUT in Prague
%
% Last update: 7/12/2016

if ~exist('bPlot', 'var')
    bPlot = 0;
end

nxData = replaceNanInFhr(nxData); 

%% modified VCH code
nxDataDtr = nxData-baseline;
fsgAccel = sgolayfilt(nxDataDtr,3,151);
fsgBeginEnd = sgolayfilt(nxDataDtr,5,151);
accPossibleParts = fsgAccel >= ampTreshold;

if bPlot
    % temporary plot
    figure(1521)
    hold off
    plot(nxDataDtr)
    hold on;
    plot(fsgBeginEnd,'r')
    x = 1:length(fsgAccel);
    %plot(x,-depthTreshold,'g','LineWidth',2)
    plot(x,fsgAccel,'--g','LineWidth',1)
    plot(x(accPossibleParts),fsgAccel(accPossibleParts),'r','LineWidth',2)
    grid on
    legend('signal','filter for begin-end','filtered for decels','possible decels') % 'treshold'
end
%accPossibleParts = (nxData-ampTreshold) >= baseline;
%accPossibleParts(nxData>erronousDataTreshold) = 0;

finalAccelAtEnd = [];
startingAccelAtBeginning = [];

accPossibleStart = find(diff(accPossibleParts)>0);
accPossibleEnd = find(diff(accPossibleParts)<0);
if length(accPossibleEnd) < length(accPossibleStart)
    finalAccelAtEnd = accPossibleStart(end);
    accPossibleStart(end)=[];
elseif length(accPossibleEnd) > length(accPossibleStart)
    startingAccelAtBeginning = accPossibleEnd(1);
    accPossibleEnd(1)=[];
end

if ~isempty(accPossibleStart),
    for i = 1:length(accPossibleStart),
        %First we have to check the length
        if (accPossibleEnd(i)-accPossibleStart(i)) > lengthTreshold*nFs,
            accStart(i) = accPossibleStart(i);
            accEnd(i) = accPossibleEnd(i);
        else
            accStart(i) = 0;
            accEnd(i) = 0;
        end
    end
    accStart(find(accStart==0))=[];
    accEnd(find(accEnd==0))=[];
    
    if isempty(accStart) && isempty(accEnd)
        accVector = nan(1,length(nxData));
        accPos = [];
        return
    end
    
    %%
    if ~isempty(finalAccelAtEnd),
        if (length(nxData) - finalAccelAtEnd) > lengthTreshold*nFs/2,
            %accVector(finalAccelAtEnd:end)=1;
            accStart(length(accStart)) = finalAccelAtEnd;
            accEnd(length(accEnd)) = length(nxData);
        end
    end
    if ~isempty(startingAccelAtBeginning)
        if startingAccelAtBeginning > lengthTreshold*nFs/2,%at least half has to be inside
            %accVector(1:startingAccelAtBeginning)=1;
            accStart = [1, accStart];
            accEnd = [startingAccelAtBeginning, accStart];
        end
    end
    %%
    accVector = NaN(1,length(nxData));
    for i = 1:length(accStart),
        accVector(accStart(i):accEnd(i))=1;
    end
    
    if ~isempty(accStart),
        for i = 1:length(accStart),
            accPos(i).start=accStart(i);
            accPos(i).end=accEnd(i);
        end
    else
        accPos = [];
    end
else
    accVector = nan(1,length(nxData));
    accPos = [];
end

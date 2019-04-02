function [decVector, decPos] = eFIGO_detectDec2(nxData, baseline, nFs, depthTreshold,lengthTreshold)

% JS: tiny modification - added fine tunning of decelerations

%%% START OF TEMPORARY CODE TO MAKE IT WORK WITH CTU-UHB Database
addedZerosAtEnd = [];
if ~isempty(find(isnan(nxData)==1)) || ~isempty(find(nxData==-1)),
    if ~isempty(find(nxData==-1)),
        addedZerosAtEnd = length(find(nxData==-1));
        nxData(nxData==-1)=[];
    end
    lengthOfData = length(nxData);
    timeOldData = 1:1:lengthOfData;
    oldData = nxData;
    timeNewData = timeOldData;
    timeNewData(isnan(nxData)==1)=[];
    nxData(isnan(nxData)==1)=[];
    data_interp = pchip(timeNewData,nxData,timeOldData);
    %figure; plot(data_interp,'r'); hold on; plot(oldData);
    nxData = data_interp;
end
if ~isempty(addedZerosAtEnd)
    clear('nxData');
    nxData = [data_interp,ones(1,addedZerosAtEnd)];
end
if size(nxData,1)~=1
    nxData = nxData';
end
%%% END OF TEMPORARY CODE TO MAKE IT WORK WITH CTU-UHB Database

%% modified VCH code
bPlot = false;

nxDataDtr = nxData-baseline;
fsgDecel = sgolayfilt(nxDataDtr,3,151);
fsgBeginEnd = sgolayfilt(nxDataDtr,5,151);

decPossibleParts = fsgDecel <= -depthTreshold;

if bPlot
    % temporary plot
    figure(1521)
    hold off
    plot(nxDataDtr)
    hold on;
    plot(fsgBeginEnd,'r')
    x = 1:length(fsgDecel);
    plot(x,-depthTreshold,'g','LineWidth',2)
    plot(x,fsgDecel,'--g','LineWidth',1)
    plot(x(decPossibleParts),fsgDecel(decPossibleParts),'g','LineWidth',2)
    grid on
    %legend('signal','filter for begin-end','treshold','filtered for decels','possible decels')
end

finalDecelAtEnd = [];
startingDecelAtBeginning = [];

decVector = nan(1,length(nxData));
decPos = [];

decPossibleStart = find(diff(decPossibleParts)>0);
decPossibleEnd = find(diff(decPossibleParts)<0);

if isempty(decPossibleStart) || isempty(decPossibleEnd),
    return;
end

if length(decPossibleEnd) < length(decPossibleStart),
    finalDecelAtEnd = decPossibleStart(end);
    decPossibleStart(end)=[];
elseif length(decPossibleEnd) > length(decPossibleStart),
    startingDecelAtBeginning = decPossibleEnd(1);
    decPossibleEnd(1)=[];
elseif decPossibleEnd(1)<decPossibleStart(1),
    decPossibleEnd(1)=[]; decPossibleStart(end)=[];
end

%% fine tune end and begin of decels
nSearch = 120*nFs;
for i = 1:length(decPossibleStart)
    
    % decel if fsgBeginEnd is above baseline
    pos = findBegin(fsgBeginEnd, decPossibleStart(i), nSearch);
    if isempty(pos)
        % gradient of decel
        pos = finetuneBegin(fsgBeginEnd,nxDataDtr,decPossibleStart(i), nSearch);
    end
    decPossibleStart(i) = pos;

    pos = findEnd(fsgBeginEnd, decPossibleEnd(i), nSearch);
    if isempty(pos)
        pos = finetuneEnd(fsgBeginEnd,nxDataDtr,decPossibleEnd(i), nSearch);
    end
    decPossibleEnd(i) = pos;
 
%     figure(1521)
%     plot(decPossibleStart(i),fsgBeginEnd(decPossibleStart(i)),'.r','MarkerSize',30);
%     plot(decPossibleEnd(i),fsgBeginEnd(decPossibleEnd(i)),'.k','MarkerSize',30);
%     
%     a = 5;
end

if bPlot
    figure(1521)
    stem(decPossibleStart,fsgBeginEnd(decPossibleStart),'.r','MarkerSize',30);
    stem(decPossibleEnd,fsgBeginEnd(decPossibleEnd),'.k','MarkerSize',30);
end

%% check overlapping decel (if overlap - join them)
for i = 2:length(decPossibleStart)
    if decPossibleStart(i) <= decPossibleEnd(i-1) 
        decPossibleEnd(i-1) = decPossibleEnd(i);
        decPossibleStart(i) = 0;
        decPossibleEnd(i) = 0;
    end
end

%% check the length of decels
for i = 1:length(decPossibleStart)
    
    if (decPossibleEnd(i)-decPossibleStart(i)) > lengthTreshold*nFs,
        decStart(i) = decPossibleStart(i);
        decEnd(i) = decPossibleEnd(i);
    else
        decStart(i) = 0;
        decEnd(i) = 0;
    end
end

decStart(decStart==0)=[];
decEnd(decEnd==0)=[];

if ~isempty(finalDecelAtEnd),
    if (length(nxData) - finalDecelAtEnd) > lengthTreshold*nFs/2,
        %        decVector(finalDecelAtEnd:end)=1;
        if ~isempty(decStart)
            decStart(length(decStart)) = finalDecelAtEnd;
            decEnd(length(decEnd)) = length(nxData);
        else
            decStart(1) = finalDecelAtEnd;
            decEnd(1) = length(nxData);
        end
        
    end
end
if ~isempty(startingDecelAtBeginning)
    if startingDecelAtBeginning > lengthTreshold*nFs/2,
        %        decVector(1:startingDecelAtBeginning)=1;
        decStart = [1, decStart];
        decEnd = [startingDecelAtBeginning, decEnd];
    end
end

decVector = NaN(1,length(nxData));
for i = 1:length(decStart),
    decVector(decStart(i):decEnd(i))=1;
end



if ~isempty(decStart),
    for i = 1:length(decStart),
        decPos(i).start=decStart(i);
        
        
        decPos(i).end=decEnd(i);
        [~,peak] = min(nxData(decStart(i):decEnd(i)));
        decPos(i).peak= decStart(i) + peak;
    end
else
    decPos = [];
end
%decPos(:).start = decStart;
%decPos.end = decEnd;

% 
% %%
% finalDecelAtEnd = [];
% startingDecelAtBeginning = [];
% 
% decPossibleStart = find(diff(decPossibleParts)>0);
% decPossibleEnd = find(diff(decPossibleParts)<0);
% 
% if ~isempty(decPossibleStart) && ~isempty(decPossibleEnd),
%     if length(decPossibleEnd) < length(decPossibleStart),
%         finalDecelAtEnd = decPossibleStart(end);
%         decPossibleStart(end)=[];
%     elseif length(decPossibleEnd) > length(decPossibleStart),
%         startingDecelAtBeginning = decPossibleEnd(1);
%         decPossibleEnd(1)=[];
%     elseif decPossibleEnd(1)<decPossibleStart(1),
%         decPossibleEnd(1)=[]; decPossibleStart(end)=[];
%     end
% end
% 
% if ~isempty(decPossibleStart) && ~isempty(decPossibleEnd)
%     for i = 1:length(decPossibleStart),
%         %First we have to check the length
%         if (decPossibleEnd(i)-decPossibleStart(i)) > lengthTreshold*nFs,
%             decStart(i) = decPossibleStart(i);
%             decEnd(i) = decPossibleEnd(i);
%         else
%             decStart(i) = 0;
%             decEnd(i) = 0;
%         end
%     end
%     decStart(decStart==0)=[];
%     decEnd(decEnd==0)=[];
%     
%     if ~isempty(finalDecelAtEnd),
%         if (length(nxData) - finalDecelAtEnd) > lengthTreshold*nFs/2,
%             %        decVector(finalDecelAtEnd:end)=1;
%             if ~isempty(decStart)
%                 decStart(length(decStart)) = finalDecelAtEnd;
%                 decEnd(length(decEnd)) = length(nxData);
%             else
%                 decStart(1) = finalDecelAtEnd;
%                 decEnd(1) = length(nxData);
%             end
%             
%         end
%     end
%     if ~isempty(startingDecelAtBeginning)
%         if startingDecelAtBeginning > lengthTreshold*nFs/2,
%             %        decVector(1:startingDecelAtBeginning)=1;
%             decStart = [1, decStart];
%             decEnd = [startingDecelAtBeginning, decEnd];
%         end
%     end
%     
%     decVector = NaN(1,length(nxData));
%     for i = 1:length(decStart),
%         decVector(decStart(i):decEnd(i))=1;
%     end
%     
%     if ~isempty(decStart),
%         for i = 1:length(decStart),
%             decPos(i).start=decStart(i);
%             decPos(i).end=decEnd(i);
%         end
%     else
%         decPos = [];
%     end
%     %decPos(:).start = decStart;
%     %decPos.end = decEnd;
% else
%     decVector = nan(1,length(nxData));
%     decPos = [];
% end

%% 
function pos = finetuneBegin(signal,nxDataDtr, iStart, nSearch)

% find proper begin
nB = iStart - nSearch;
if nB <= 0,
    nB = 1;
end

sig = signal(nB:iStart);
data = nxDataDtr(nB:iStart);
sig2 = slidingAvg(sig,50);

% figure(1)
% hold off
% plot(sig)
% hold on;
% plot(sig2,'r');
% title('BEGIN')

for t = length(sig2):-1:2
    if data(t) > 0 || sig2(t-1) < sig2(t)
        break;
    end
end

% plot(t,sig2(t),'.r','MarkerSize',30);
if isempty(t)
    t=0;
end

pos = nB + t;

function pos = findBegin(signal, iStart, nSearch)

% find proper begin
nB = iStart - nSearch;
if nB <= 0,
    nB = 1;
end

idx = find(signal(nB:iStart) >= 0);
% 
% figure(1)
% hold off
% plot(signal(nB:iStart),'b')
% title('BEGIN')

if ~isempty(idx)
    pos = nB + idx(end);
else
    pos = [];
end

%%
function pos = finetuneEnd(signal,nxDataDtr, iEnd, nSearch)

% find proper begin
nE = iEnd + nSearch;
if nE > length(signal)
    nE = length(signal);
end

sig = signal(iEnd:nE);
data = nxDataDtr(iEnd:nE);
sig2 = slidingAvg(sig,50);

p2 = polyfit(1:length(sig2),sig2,1);
bestFit = polyval(p2,1:length(sig2));

% figure(2)
% hold off
% plot(data,'b')
% hold on;
% plot(sig,'r')
% plot(sig2,'m');
% plot(bestFit,'k')
% title('END')

for t = 1:length(sig2)-1
    if data(t) > 0 || sig2(t) > sig2(t+1)
        break;
    end
end

% plot(t,sig2(t),'.k','MarkerSize',30);

pos = iEnd + t;

%%
function pos = findEnd(signal, iEnd, nSearch)

nE = iEnd + nSearch;
if nE > length(signal)
    nE = length(signal);
end

idx = find(signal(iEnd:nE) >= 0);

% figure(2)
% hold off
% plot(signal(iEnd:nE),'b')
% title('END')

if ~isempty(idx)
    pos = iEnd + idx(1) - 1;
else
    pos = [];
end

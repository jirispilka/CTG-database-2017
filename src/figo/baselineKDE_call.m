function [baselineOut,timePointS,baselinePointS,baselinePointS_unsure, baselinePointS_unassigned] = baselineKDE_call(nxData,nFs,nSegmentLength, nOverlap, bPlot)
% [baselineOut,timePointS,baselinePointS,baselinePointS_unsure, baselinePointS_unassigned] = baselineKDE_call(nxData,nFs,nSegmentLength, nOverlap, bPlot)
%
% Baseline computed on basis of paper by Georgieva 2011. The method uses
% segments of signal. Default segment value is 5 minutes with 50% overlap.
%
% INPUTS:
%       nxData          -  [bpm] data vector
%       nFs             -  [Hz] sampling frequency of nxData
%       nSegmentLength  -  [min] length of segment in which to compute the baseline
%                          (default 5 minutes)
%       nOverlap        -  [%] window overlap
%       bPlot           -  [boolean] should the function plot the outcome
%
% OUTPUTS:
%       baselineOut             - baseline fitted with cubicspline
%       timePointS              - position of baseline points in time
%       baselinePointS          - points of baseline
%       baselinePointS_unsure   - points of "not sure" baseline
%       baselinePointS_unassigned   - points of "unasignable" baseline
%
% Authors:
%       Vaclav Chudacek, Lukas Zach
%
% Last update: 19/12/2016
%
% Reference:
%       Title:      Computerized fetal heart rate analysis in labor:
%                   detection of intervals with un-assignable baseline
%       Authors:    Antoniya Georgieva et al.
%       BibName:    Georgieva2011

%% In case only stable segments should be used set useStableSegments = 1;
useOnlyStableSegments = 0;

%% Input arguments revised
if nargin < 2
    error('Sampling frequency of the data is missing!');
elseif nargin < 3
    nSegmentLength = 5;
elseif nargin < 4
    %error('You have to put in either both segment length and overlap, or none of them!');
    nOverlap = 50; % Default value 50%
elseif nargin < 5
    bPlot = 0;
end
nSegmentLength = nFs * 60 * nSegmentLength; % Default value 300s

if length(nxData) < nSegmentLength
    error('The data provided are too short for baseline estimation. Change the length of the segment, or input longer data-vector.');
end
if nOverlap~=50 && nOverlap~=75
    disp('Only 50/75% overlap is supported. Changing to default 50%.')
    nOverlap = 50;
end

%%
% Constants
overlapRatio = 100/(100-nOverlap);
partialWindowLength = nSegmentLength/(overlapRatio);
numOfWindows = fix(length(nxData)/partialWindowLength);

for j = 1 : numOfWindows % Preparation of sub-windows for computation
    actualWindow = split2segments(j,overlapRatio,numOfWindows,nxData,partialWindowLength,nSegmentLength);
    lengthOfData = length(actualWindow);
    indxMiss1 = find(actualWindow==-1);
    indxMiss2 = find(isnan(actualWindow)==1);
    indxMiss3 = find(actualWindow==0); 
    %indxMiss3 = find(actualWindow<50); % alternative selection of data with erronous values 
    missRatio = (length(indxMiss1)+length(indxMiss2)+length(indxMiss3))/lengthOfData;
    if missRatio > 0.49
       confidenceInData(j) = -1;
    else
        confidenceInData(j) = 1;
    end
end

%% We need to prepare data when nans or -1 are present 
% % Additional pre-processing ENSL_db - outcome is used at the end of the file
nanCut = 0;
if isnan(nxData(end-100:end)),
    indexLastNaN = find(isnan(nxData(end-100:end))==1,1,'last'); 
    nxData(end)=nanmedian(nxData);
end

% (Pre-processing for the CTU-UHB)
addedZerosAtEnd = [];
if ~isempty(isnan(nxData)==1) || ~isempty(find(nxData==-1,1)),
    disp('Preprocessing... -> Interpolating missing values.')
    if ~isempty(find(nxData==-1,1)),
        addedZerosAtEnd = length(find(nxData==-1));
        nxData(find(nxData==-1,1))=[];
    end
    lengthOfData = length(nxData);
    timeOldData = 1:1:lengthOfData;
    %oldData = nxData;
    timeNewData = timeOldData;
    timeNewData(isnan(nxData)==1)=[];
    nxData(isnan(nxData)==1)=[];
    % When nans are present in the signal
    timeOldData_max = find(timeOldData > max(timeNewData),1);
    if isempty(timeOldData_max)
        timeOldData_max = length(timeOldData);
    end
    %nanCutNmbr = length(timeOldData)-timeOldData_max;
    % We resample the data without the nans
    data_interp = pchip(timeNewData,nxData,timeOldData(1:timeOldData_max));
    clear('nxData');
    nxData = data_interp;
    disp(['Computing... -> Baseline in window of ',num2str(nSegmentLength/(nFs * 60)),' minutes.'])
end
% if bPlot,
%     figure; plot(oldData); hold on; plot(data_interp,'r'); 
% end

%% Computing baseline points
% Set up of variables
baselinePointS = zeros(numOfWindows,1);
baselinePointS_unsure = zeros(numOfWindows,1);
baselinePointS_unassigned = zeros(numOfWindows,1);
timePointS = zeros(numOfWindows,1);
timePointS_unnasigned = zeros(numOfWindows,1);


%% Segmentation 

for j = 1 : numOfWindows % Preparation of sub-windows for computation
   actualWindow = split2segments(j,overlapRatio,numOfWindows,nxData,partialWindowLength,nSegmentLength);
    % We take out all zeros for computation of the baseline
    actualWindow(actualWindow==0)=[];
    % If full zeros or constant...
    if isempty(actualWindow) || (min(actualWindow)-max(actualWindow) == 0)
        resultStability = 4; % No reliable information at all - either empty or precomputed...
        if isempty(actualWindow)
            baselinePoint = -1;
        else
            baselinePoint = actualWindow(1);
        end
    else
        % Running the main file computing the KDE and find stable seg.
        n=128;
        aprioriConf = confidenceInData(j);
        [SSI,resultStability,baselinePoint,bandwidth,density,xmesh] = baselineKDE_comp(actualWindow,aprioriConf,n);
%         if bPlot,
%             figure(2000);
%             hold on;
%             plot(xmesh,density,'r');
%             figure(3000);
%             hold on;
%             plot(xmesh(1:length(fastmovingaverage(density,5))),fastmovingaverage(density,5),'k');
%         end
        % NOT TESTED!!! Additional information about stable segments in the sub-window
        if useOnlyStableSegments
            stableSegm = findStableSegment(actualWindow,nFs,0);
            baselinePointSS(j) = median(actualWindow(stableSegm));
            if bPlot && 1==0
                plot(1:2000,baselinePoint,'g');
                plot(1:2000,median(actualWindow(stableSegm)),'r');
            end
        else
            baselinePointSS = [];
        end
    end
    %Division of baseline points according to assumed stability of the baseline
    if resultStability==1
        baselinePointS(j) = baselinePoint;
        timePointS(j) = j*partialWindowLength;
    elseif resultStability==2
        baselinePointS(j) = baselinePoint;
        timePointS(j) = j*partialWindowLength;
    else
        baselinePointS_unassigned(j) = baselinePoint;
        if ~isempty(actualWindow),
            timePointS_unnasigned(j) = j*partialWindowLength;
        end
    end
end

%% Cleaning up the points
% We have to clear unnasignable baseline points from the vectors
genTimePoint = timePointS + timePointS_unnasigned;
timePointS(find(baselinePointS==0))=[];
baselinePointS(find(baselinePointS==0))=[];

if isempty(baselinePointS) % In case we have no hit - we use median
%     baselineOut(1:length(nxData)) = median(nxData);
    finalBaseline(1:length(nxData)) = median(nxData);
else
    % we add first point for smoother interpolation
    baselinePointS = [baselinePointS(1); baselinePointS];
    timePointS = [1; timePointS];
    % we add the last point if not there already
    timePointS = [timePointS; length(nxData)];
    if length(timePointS) == length(unique(timePointS))
        baselinePointS = [baselinePointS; baselinePointS(end)];
    else
        timePointS = unique(timePointS);
    end
    % we interpolate to get the final baseline
    % pchip seems to work better... than spline
    finalBaseline = pchip(timePointS,baselinePointS,1:1:timePointS(end));
    if ~isempty(addedZerosAtEnd)
        finalBaseline = [finalBaseline,ones(1,addedZerosAtEnd)];
    end
end
baselineOut = finalBaseline;

if bPlot
    if ~exist('baselinePointSS', 'var')
        baselinePointSS = [];
    end
    plotBaselineAccDec(nxData,baselineOut,nFs,genTimePoint,timePointS,baselinePointS,baselinePointSS,'k');
end

% if nanCut,
%     baselineOut =  [finalBaseline,nan(1,nanCutted)];
% end

function actualWindow = split2segments(j,overlapRatio,numOfWindows,nxData,partialWindowLength,nSegmentLength)
if overlapRatio > 1
    ind = (j-1) * partialWindowLength;
    if j > numOfWindows - (overlapRatio - 1)
        actualWindow = nxData(ind+1:end);
    elseif ind + nSegmentLength > length(nxData)
        actualWindow = nxData(ind+1:end);
    else
        actualWindow = nxData(ind+1:ind + nSegmentLength);
    end
else
    %NOT TESTED! No-overlap in the data
    actualWindow = nxData((j-1)*partialWindowLength+1:(j-1)*partialWindowLength+nSegmentLength);
end

function [cFeatures] = uti_computeFigoSegment(aData,aBegEnd,eFIGO,fs)

bMakeAllNaN = false;
if isnan(aData)
    % compute with random and make all NaN at the end
    aData = randn(size(aData));
    bMakeAllNaN = true;
end

% [out] = analyzeFHR_enhancedFIGO(aData',fs,'',0);

nBegin = aBegEnd(:,1);
nEnd = aBegEnd(:,2);
eFIGO.baseLine = eFIGO.baseLine(nBegin:nEnd);

decel_vec = eFIGO.decels.vector(nBegin:nEnd);
xaxis = 1:length(decel_vec);

% figure
% plot(aData)
% hold on;
% plot(eFIGO.baseLine,'m')
% plot(xaxis(decel_vec==1),aData(decel_vec==1),'r')
% %plot(eFIGO.decels(nBegin:nEnd),'m')

cFeatures.bslnMean = nanmean(eFIGO.baseLine);
cFeatures.bslnSD = nanstd(eFIGO.baseLine);
cFeatures.stressRatio = nansum(decel_vec)/length(decel_vec);

% accels
if ~isempty(eFIGO.accels.position)
    acc_b = [eFIGO.accels.position.start];
    acc_e = [eFIGO.accels.position.end];
    cFeatures.accNumber = sum(acc_b > nBegin & acc_e < nEnd);
else
    cFeatures.accNumber = 0;
end

% decels
dec_segment = [];
if ~isempty(eFIGO.decels.position)
    dec_pos = eFIGO.decels.position;
    cnt = 0;
    for i = 1:length([dec_pos.start])
        if dec_pos(i).start >= nBegin && dec_pos(i).end < nEnd
            cnt = cnt + 1;
            dec_segment(cnt).start = dec_pos(i).start - nBegin + 1;
            dec_segment(cnt).end = dec_pos(i).end - nBegin + 1;
            dec_segment(cnt).peak = dec_pos(i).peak - nBegin + 1;
        end
    end
    
    if isempty(dec_segment)
        cFeatures.decNumber = 0;
        cFeatures.areaDecelTriangle = 0;
    else
        cFeatures.decNumber = length([dec_segment.start]);
        area_tri_samp = eFIGO_analyzeDecels(aData,dec_segment);
        cFeatures.areaDecelTriangle = area_tri_samp;
    end
else
    cFeatures.decNumber = 0;
    cFeatures.areaDecelTriangle = 0;    
end

% STV, LTV
cFeatures.LTV_FIGO_bpm = featureLTV_FIGO(aData.*eFIGO.multVect(nBegin:nEnd)',fs);
cFeatures.STV_Sonicaid = featureSTV_Sonicaid(aData.*eFIGO.multVect(nBegin:nEnd)',fs);

% computations
idxNotNan = ~isnan(eFIGO.baseLine);
xx = 1:length(eFIGO.baseLine);
p = polyfit(xx(idxNotNan),eFIGO.baseLine(idxNotNan),1);

cFeatures.bslnBeta0 = p(2);
cFeatures.bslnBeta1 = p(1);

if ~isempty(dec_segment)
    x = [dec_segment.peak];
    x = diff(x)/fs;
    if isempty(x)
        x = 0;
    end
else
    x = 0;
end

cFeatures.decDeltaMedian = nanmedian(x);
cFeatures.decDeltaMad = mad(x);

dtrd = aData - eFIGO.baseLine';
cFeatures.decDtrdPlus = sum(dtrd(dtrd > 0))/length(dtrd);
cFeatures.decDtrdMinus = sum(dtrd(dtrd < 0))/length(dtrd);
cFeatures.decDtrdMedian = nanmedian(dtrd);
cFeatures.decDtrdMAD = mad(dtrd);

%% MAKE ALL FEATURES NAN
if bMakeAllNaN
    cnames = fieldnames(cFeatures);
    for i = 1:length(cnames)
        cFeatures.(cnames{i}) = NaN;
    end
end

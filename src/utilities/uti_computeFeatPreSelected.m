function [cFeatures] = uti_computeFeatPreSelected(aData,sSelFeat,aBegEnd,fs,eFigo)

cFeatures = [];

bMakeAllNaN = false;
if isnan(aData)
    % compute with random and make all NaN at the end
    aData = randn(size(aData));
    bMakeAllNaN = true;
end

aDataPossibleNaN = aData;
aData = removeNaNsAtBeginAndEnd(aData);

sSelFeat = lower(sSelFeat);

%% statistics
if contains(sSelFeat, 'stat')
    cFeatures.mean = nanmean(aData);
    cFeatures.std = nanstd(aData);
    cFeatures.median = nanmedian(aData);
    idx = ~isnan(aData);
    cFeatures.mad = mad(aData(idx));
    cFeatures.skewness = skewness(aData(idx));
    cFeatures.kurtosis = kurtosis(aData(idx));
    %cFeatures.quality = 1-sum(isnan(aData))/length(aData);
end

%% FIGO
if contains(sSelFeat, 'figo') || contains(sSelFeat,'all')
    
    featfigo = uti_computeFigoSegment(aDataPossibleNaN,aBegEnd,eFigo,fs);
    
    % the trend is intersting on the whole record not on the segment
    featfigo.bslnAllBeta0 = eFigo.bslnAllBeta0;
    featfigo.bslnAllBeta1 = eFigo.bslnAllBeta1;
    
    cFeatures = copyStruct(featfigo, cFeatures);
end

%% spectral
if contains(sSelFeat,'spectral') || contains(sSelFeat,'all')
    
    % interpolate all
    tempData = interpolateAllGaps(aData, fs, 0);
    [~,~,~,~,energySpec,alpha] = HFLFanalysisSplit(tempData,fs,[0.04 1],[0.04 1],[], 0);
        
    % based on Siira 2013
    cFeatures.energy_VLF = energySpec(1)+energySpec(2);
    cFeatures.energy_LF = energySpec(3);
    cFeatures.energy_HF = energySpec(4) + energySpec(5);
    cFeatures.energy_LF_HF = cFeatures.energy_LF/cFeatures.energy_HF;
    cFeatures.energy_tot = energySpec(6);
    cFeatures.spectrum_slope = alpha;
end


%% MF and H
if contains(sSelFeat,'mfhurst') || contains(sSelFeat,'all')
      
    gamint = 1;
    j1 = 3;
    j2 = 10;
    %cfeat1 = featuresMF_LeadersMulti_simple(aData',fs,gamint);
    [cfeat,~,logstat] = featuresMF_RFL(aData',j1,j2,gamint);
    cFeatures = copyStruct(cfeat,cFeatures);
    cFeatures.H310 = cfeat.MF_HDWT;
    cFeatures = rmfield(cFeatures,'MF_HDWT');
       
    j1 = 2;
    j2 = 9;
    cfeat = featuresMF_RFL(logstat,j1,j2,gamint);
    cFeatures.H29 = cfeat.MF_HDWT;
    cFeatures.MF_c1_29 = cfeat.MF_c1;
    cFeatures.MF_c2_29 = cfeat.MF_c2;
    cFeatures.MF_c3_29 = cfeat.MF_c3;
    cFeatures.MF_c4_29 = cfeat.MF_c4;
    
end

%% MAKE ALL FEATURES NAN
if bMakeAllNaN
    cnames = fieldnames(cFeatures);
    for i = 1:length(cnames)
        cFeatures.(cnames{i}) = NaN;
    end
end

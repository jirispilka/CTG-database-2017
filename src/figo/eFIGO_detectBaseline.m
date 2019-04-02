function [newBl_diag, newBl_accDec, basal] = eFIGO_detectBaseline(fhr,fs,bPlotAlgoSteps)

% [diagBaaseline_out, decelBaseline_out] = eFIGO_detectBaseline(fhr,fs,bPlotAlgoSteps)
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
%       Vaclav Chudacek, Jiri Spilka
%       http://ctg.ciirc.cvut.cz
%       CIIRC, CVUT in Prague
%
% Last update: 19/12/2016

%% Signal is undersampled on places with potentially deep decelerations (based on a threshold)
underSampAmpLimit = 100; % upperlimit of FHR where undersampling will be performed
underSampFactor = 2; % undersampling factor to be used for signal culling

fhr_underSampled = fhr;
indx_temp = find(fhr_underSampled < underSampAmpLimit);
indx = indx_temp(1:underSampFactor:end);
fhr_underSampled(indx)=[];

%% Estimation of BASELINEs with different window length based on undersampled FHR

bl1 = baselineKDE_call(fhr_underSampled, fs, 20   , 50, bPlotAlgoSteps);
bl2 = baselineKDE_call(fhr_underSampled, fs,  5   , 50, bPlotAlgoSteps);
bl3 = baselineKDE_call(fhr_underSampled, fs,  2.5 , 50, bPlotAlgoSteps);
bl4 = baselineKDE_call(fhr_underSampled, fs,  1.25, 50, bPlotAlgoSteps);

if bPlotAlgoSteps
    figure(1000); clf; grid on; hold on;
    plot(fhr);
end

bl_accDec = max([bl1;bl2;bl3]); % Estimation of accDec baseline is based on max of different baselines [bl1 - bl3]
bl_diag = fastmovav(bl_accDec, fs*20*60+1); % smoothed 20-min window baseline is used as a diagnostic bl

if bPlotAlgoSteps
    figure(1000); clf; grid on; hold on;
end

%% Estimation of (correct-length) baselines
newBl_accDec = bl_accDec;
for i = 1:length(indx)
    indxNew(i) = indx(i)+i-1;
    newBl_accDec = [newBl_accDec(1:min(length(newBl_accDec)-1,indxNew(i))),newBl_accDec(min(length(newBl_accDec)-1,indxNew(i))),newBl_accDec(min(length(newBl_accDec)-1,indxNew(i))+1:end)];
end
newBl_diag = fastmovav(newBl_accDec, fs*20*60+1);

% N=fs*20*60;
% newBl_diag = filter(ones(N,1)/N,1,newBl_accDec);

if bPlotAlgoSteps
    figNum = 1000; length2plot = 10000;
    plotBaseline(figNum, length2plot,newBl_accDec,newBl_diag,bl1,bl2,bl3,bl4,fhr,fs)
    figNum = 2000;  length2plot = length(newBl_accDec);
    plotBaseline(figNum, length2plot,newBl_accDec,newBl_diag,bl1,bl2,bl3,bl4,fhr,fs)
end

%% BASAL line
basal = zeros(length(fhr),1);
samp_min = 60*fs;
% baseline based approach
confBasal = 'low';
basal(:,1) = nanmedian(newBl_diag(1*samp_min:10*samp_min)',1)';
if basal(1,1) < 100 || basal(1,1) > 165
    basal(:,1) = nanmedian(newBl_diag(10*samp_min:20*samp_min)',1)';
end

function plotBaseline(figNum, length2plot,bl_accDec,bl_diag,bl1,bl2,bl3,bl4,fhr,fs)
%% Plotting
figure(figNum); hold on;
timeVect = 1/60/fs:1/60/fs:length(bl_accDec(1:length2plot))/60/fs;
plot(timeVect,fhr(1:length2plot));
h_accDec = plot(timeVect,bl_accDec(1:length2plot),'y','LineWidth',5);
h_diag = plot(timeVect,bl_diag(1:length2plot),'--m','LineWidth',5);
h_diagBe = plot(timeVect,bl_diag(1:length2plot)-5,'-m','LineWidth',2);
h_diagAb = plot(timeVect,bl_diag(1:length2plot)+5,'-m','LineWidth',2);
plot(timeVect(1:min(length(bl1),length2plot)),bl1(1:min(length(bl1),length2plot)),'m'); 
plot(timeVect(1:min(length(bl1),length2plot)),bl2(1:min(length(bl2),length2plot)),'g'); 
plot(timeVect(1:min(length(bl1),length2plot)),bl3(1:min(length(bl3),length2plot)),'k'); 
plot(timeVect(1:min(length(bl1),length2plot)),bl4(1:min(length(bl4),length2plot)),'y');

legend('FHR','blAccDec','blDiag', 'upperLimit', 'lowerLimit','20-min', '5-min', '2.5-min', '1.25-min');




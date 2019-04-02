function [] = plotBaselineAccDec(nxData,baselineOut,nFs,genTimePoint,timePointS,baselinePointS,baselinePointSS,baseColor)
figure(1000); hold on;
timeVect = 1/60/nFs:1/60/nFs:length(nxData)/60/nFs;
tempData = nxData;
nxData(nxData<30)=nan;
h0 = plot(timeVect,nxData); hold on;
h1 = plot(timeVect,baselineOut,baseColor,'LineWidth',5);

plot(timePointS/60/nFs,baselinePointS,'xr');
genTimePoint(genTimePoint==0)=[];
if ~isempty(baselinePointSS),
    plot(genTimePoint/60/nFs,baselinePointSS(1:length(genTimePoint)),'ok');
    legend('FHR','Computed baseline','EST:Kernel','EST:Stable segm.','Location','SouthWest')
else
    legend('FHR','Computed baseline','EST:Kernel','Location','SouthWest')
end

nxData = tempData;
xlabel('Time [min]');
ylabel('FHR [bpm]');
title('Example of FHR-baseline estimation on recording from CTU-UHB db')
%         [decVec,accVec] = detectAccelDecel(nxData,4,'georgieva');
%         figure(1000);
%         h2 = plot(timeVect,decVec*100,'xr');
%         plot(timeVect,accVec*150,'xg');
%         legend([h0,h1,h2],'FHR','Estimated baseline','Data-points belonging to deceleration','Location','SouthWest')

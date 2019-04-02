function [area_tri_samp2,decel] = eFIGO_analyzeDecels(fhr,decel)
%%
% stressRatio = 0;
area_tri_samp2 = 0;

if isempty(decel)
    return
end

%%
% stressTime = 0;
% for i = 1:length(decel),
%     stressTime = stressTime + (decel(i).end - decel(i).start);
% end
% 
% %nonStressTime = length(fhr)-stressTime;
% %stressRatio = nonStressTime/stressTime;
% 
% stressRatio = stressTime/length(fhr);

%% analyse decel by area (triangular)

%area_tri_samp2 = zeros(length(decel),1);

for i = 1:length(decel),
    nB = decel(i).start;
    nE = decel(i).end;
    
    baseline = nanmean([fhr(nB),fhr(nE)]);
    [nMinDecel, nMinInd] = min(fhr(nB:nE));
    nAmp = abs(baseline - nMinDecel);
    
    decel(i).area_triangle = ((nE - nB)*nAmp)/2;
end

area_tri_samp2 = nanmedian([decel.area_triangle]);

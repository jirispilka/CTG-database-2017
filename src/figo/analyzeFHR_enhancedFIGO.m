function [out] = analyzeFHR_enhancedFIGO(fhr,fs,fhrID,verbose)

%% General parameters
unitShiftDiff = 5; % setup of the unit [bpm] for "confidence lines" of baseline 

if verbose > 1
    bPlot = 1;
    bPlotAlgoSteps = 1; % all algorithm steps will be plotted out
elseif verbose == 1
    bPlot = 1;
    bPlotAlgoSteps = 0; 
else
    bPlot = 0;
    bPlotAlgoSteps = 0; 
end

%%
 [newBl_diag, newBl_accDec] = eFIGO_detectBaseline(fhr,fs,bPlotAlgoSteps);

%%
% Detection of shifts in baseline (together with "confidence" lines)
[correctedNewBl_accDec,shiftPos,shiftDirs,shiftStartPoint,shiftEndPoint,shifSigDiff] = eFIGO_detectBaseLineShifts(fhr, newBl_diag, newBl_accDec,fs,unitShiftDiff);

% %%
% %TEMP CODE
% fhr_binary = fhr;
% fhr_binary(fhr_binary<61) = 0;
% fhr_binary(fhr_binary>60) = 1;
% template = strel('line', 20, 0);
% 
% fhr_binary = imerode(fhr_binary,template);
% 
% fhr_binary = bwareaopen(fhr_binary,50);
% fhr(fhr_binary==0)=0;
% 
% 
% 
% 
% zeroLength = 4;
% gapDist = 80;
% end0 = find(fhr'~=0,1,'last');
% fhr1 = fhr;
% fhr1(end0:end)=[];
% vec_out = strfind(fhr1', zeros(1,zeroLength));
% for i=length(vec_out):-1:2,
%     if vec_out(i)-1 == vec_out(i-1),
%         vec_out(i)=[];
%     end;
% end
% 
% fhr2 = fliplr(fhr1');
% vec_out2 = strfind(fhr2, zeros(1,zeroLength));
% for i=length(vec_out2):-1:2,
%     if vec_out2(i)-1 == vec_out2(i-1),
%         vec_out2(i)=[];
%     end;
% end
% vec_out2 = abs(vec_out2 - length(fhr2));
% vec_out2 = fliplr(vec_out2);
% %diff(vec_out2)
% 
% figure;
% fhr_temp = fhr;
% fhr_temp(fhr_temp<60)=NaN;
% plot(fhr_temp);
% hold on;
% 
% vec_out(find(diff(vec_out)<gapDist)+1)=[];
% vec_out2(diff(vec_out2)<gapDist)=[];
% 
% counter=0;delete=[];
% for i = 1:length(vec_out)-1,
%     closest2 = find(vec_out2>vec_out(i),1,'first');
%     if ~isempty(closest2)
%         if vec_out2(closest2)-vec_out(i)>vec_out(i+1)-vec_out(i)
%             counter=counter+1;
%             delete(counter) = i+1;
%         end
%     end
% end
% vec_out(delete)=[];
% 
% counter=0;delete2=[];
% for i = 1:length(vec_out2)-1,
%     closest2 = find(vec_out>vec_out2(i),1,'first');
%     if ~isempty(closest2)
%     if vec_out(closest2)-vec_out2(i)>vec_out2(i+1)-vec_out2(i)
%         counter=counter+1;
%         delete2(counter) = i;
%     end
%     end
% end
% vec_out2(delete2)=[];
% 
% 
% for i = 1:length(vec_out),
%     plot(vec_out(i),80:170,'go');
% end
% 
% for i = 1:length(vec_out),
%     plot(vec_out(i),120:170,'yo');
% end
% 
% for i = 1:length(vec_out2),
%     plot(vec_out2(i),80:170,'ro');
% end
% 
% 
% for i = 1:length(vec_out2),
%     plot(vec_out2(i),120:170,'ko');
% end
% 
% figure;
% plot(fhr); hold on;
% fhr_binary = fhr;
% fhr_binary(fhr_binary<61) = 0;
% fhr_binary(fhr_binary>60) = 1;
% plot(fhr_binary*50);
% template = strel('line', 20, 0);
% fhr_binary = imerode(fhr_binary,template);
% %fhr(fhr_binary==0)=0;
% plot(fhr_binary*40,'k');
% trt = bwareaopen(fhr, 20);
% plot(trt,'g')
% 
% for i = 1:length(vec_out2),
%     fhr(vec_out(i):vec_out2(i)+1) = NaN;
% end
% plot(fhr,'r');

%%
%YOU CAN DO SPLIT TO THE SUBWINDOWS HERE!!!

%% Decelerations / Accelerations
depthTreshold = 15; lengthTreshold = 10;
%[decVec,decPos] = eFIGO_detectDec(fhr, correctedNewBl_accDec, fs, depthTreshold,lengthTreshold);
[decVec,decPos] = eFIGO_detectDec2(fhr, correctedNewBl_accDec, fs, depthTreshold,lengthTreshold);
depthTreshold = 15; lengthTreshold = 15;
[accVec,accPos] = eFIGO_detectAcc(fhr, correctedNewBl_accDec, fs, depthTreshold,lengthTreshold);

% Additional parameters from analysis of Decelerations/Accelerations
[bradyVec, bradyPos] = eFIGO_detectBrady(fhr,fs,decPos);
stressRatio = nansum(decVec)/length(decVec);
[area_tri_samp2,decPos] = eFIGO_analyzeDecels(fhr,decPos);
%reactivity = eFIGO_detectReactivity(fhr,fs,accPos);

%% STV / LTV / Baseline values
multVect = ones(1,length(decVec)); multVect(decVec==1)=nan; multVect(accVec==1)=nan;% we get rid of all parts of signal with Dec/Acc for computation of STV and LTV
if size(multVect,1)>1
    multVect = multVect';
end
[STV_Sonicaid] = featureSTV_Sonicaid(fhr.*multVect, fs, 'US');
[STV_Sonicaid_withDecel] = featureSTV_Sonicaid(fhr, fs, 'US');
[LTV_FIGO] = featureLTV_FIGO(fhr.*multVect, fs);
[LTV_FIGO_withDecel] = featureLTV_FIGO(fhr, fs);
medianBaseLine = median(newBl_diag);

%% Just plot of the whole baseline estimation process
if bPlot,
    figure(1001); clf; hold on; grid on; %set(gca,'FontSize',28);
    timeVect = 1/60/fs:1/60/fs:length(correctedNewBl_accDec)/60/fs;
    fhr(fhr==0)=nan;
    h_fhr = plot(timeVect,fhr);
    h_diag = plot(timeVect,newBl_diag,'--m','LineWidth',5);
    
    if ~isempty(shiftPos),
        for i = 1:length(shiftPos),
            if shiftDirs(i)==1,
                hshUp = plot(timeVect(shiftPos(i)),140:5:180,'^r','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');
            else
                hshDo = plot(timeVect(shiftPos(i)),120:5:160,'vb','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b');
            end
        end
    end
    idx = ~isnan(decVec);
    h_diagConf = plot(timeVect,newBl_diag-shifSigDiff*unitShiftDiff,'-.m','LineWidth',2);
    h_diagConf = plot(timeVect,newBl_diag+shifSigDiff*unitShiftDiff,'-.m','LineWidth',2);
    h_accDec = plot(timeVect,correctedNewBl_accDec,'y','LineWidth',5);
    h_dec = plot(timeVect,decVec*100,'xr','MarkerSize',8);
    h_dec = plot(timeVect(idx),fhr(idx),'r');
    h_acc = plot(timeVect,accVec*180,'xg','MarkerSize',8);
    h_brady = plot(timeVect,bradyVec*99,'+k','MarkerSize',2); plot(timeVect,bradyVec*101,'+k','MarkerSize',2);
    title([fhrID,' [StressRatio = ',num2str(stressRatio,'%1.2f'),'] with Baseline = ',num2str(round(medianBaseLine)),' bpm. | STV = ',num2str(round(STV_Sonicaid)),'ms (', num2str(round(STV_Sonicaid_withDecel)),'ms) | LTV = ',num2str(round(LTV_FIGO)),'bpm (',num2str(round(LTV_FIGO_withDecel)),'bpm)']);
    legend([h_accDec,h_diag,h_diagConf,h_acc,h_dec,h_brady],'Baseline used for Acc/Dec computation','Baseline used for Diagnostics purposes','"Confidence" limit for baseline','Accelerations','Decelerations','Prolonged Bradycardia','Location','Best');
end

%% We formate the output variable (TEMPORARY, should reflect windows)
out.baseLine = newBl_diag;
out.baseLineAccDecc = correctedNewBl_accDec;
out.decels.vector = decVec;
out.decels.position = decPos;
out.accels.vector = accVec;
out.accels.position = accPos;
out.brady.vector = bradyVec;
out.brady.position = bradyPos;
out.stressRatio = stressRatio;
out.areaDecelTriangle = area_tri_samp2;
out.STV_Sonicaid = STV_Sonicaid;
out.LTV_FIGO = LTV_FIGO;
out.multVect = multVect;

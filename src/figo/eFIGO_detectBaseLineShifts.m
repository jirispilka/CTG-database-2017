function [bl_accDec,shiftPos,shiftDirs,shiftStartPoint,shiftEndPoint,shifSigDiff] = eFIGO_detectBaseLineShifts(fhr, bl_diag, bl_accDec,fs,unitShiftDiff)

segmentLength = 5;
%shiftUpOut = []; shiftDownOut= [];
noShift = 0;
shiftPos = []; shiftDirs = [];
lastShiftPos = 1; lastShiftDir = 0; lastTriedPos = [];
meanBL = mean(bl_diag(1:segmentLength*60*fs));
maxSearchLength = 10*60*fs;
maxSearchLength2 = 2*maxSearchLength;

%%
changeUp = find(bl_diag>1.1*meanBL,1,'first');
changeDown = find(bl_diag<0.9*meanBL,1,'first');
while noShift == 0,
    if isempty(lastShiftPos)
        changeUp = find(bl_diag(lastTriedPos:end)>1.1*meanBL,1,'first')+lastTriedPos;%-1;
        changeDown = find(bl_diag(lastTriedPos:end)<0.9*meanBL,1,'first')+lastTriedPos;%-1;
    else
        changeUp = find(bl_diag(lastShiftPos:end)>1.1*meanBL,1,'first')+lastShiftPos;%-1;
        changeDown = find(bl_diag(lastShiftPos:end)<0.9*meanBL,1,'first')+lastShiftPos;%-1;
    end
    if isempty(changeDown),
        if isempty(changeUp),
            noShift = 1;
        else
            [lastShiftPos, lastShiftDir, meanBL2] = checkUp(bl_diag, meanBL, changeUp, segmentLength, fs);
            if isempty(lastShiftPos)
                if isempty(lastTriedPos),
                    lastTriedPos = changeUp;
                else
                    if lastTriedPos + fs*10 > changeUp,
                        lastTriedPos = changeUp;
                    else
                        lastTriedPos = changeUp + fs*15;
                    end
                end
            end
        end
    else
        if isempty(changeUp),
            [lastShiftPos, lastShiftDir, meanBL2] = checkDown(bl_diag, meanBL, changeDown, segmentLength, fs);
            if isempty(lastShiftPos)
                if isempty(lastTriedPos),
                    lastTriedPos = changeDown;
                else
                    if lastTriedPos + fs*10 > changeDown,
                        lastTriedPos = changeDown;
                    else
                        lastTriedPos = changeDown + fs*15;
                    end
                end
            end
        else
            if changeUp < changeDown
                [lastShiftPos, lastShiftDir, meanBL2] = checkUp(bl_diag, meanBL, changeUp, segmentLength, fs);
                if isempty(lastShiftPos)
                    if isempty(lastTriedPos),
                        lastTriedPos = changeUp;
                    else
                        if lastTriedPos + fs*10 > changeUp,
                            lastTriedPos = changeUp;
                        else
                            lastTriedPos = changeUp + fs*15;
                        end
                    end
                end
            else
                [lastShiftPos, lastShiftDir, meanBL2] = checkDown(bl_diag, meanBL, changeDown, segmentLength, fs);
                if isempty(lastShiftPos)
                    if isempty(lastTriedPos),
                        lastTriedPos = changeDown;
                    else
                        if lastTriedPos + fs*10 > changeDown,
                            lastTriedPos = changeDown;
                        else
                            lastTriedPos = changeDown + fs*15;
                        end
                    end
                end
            end
        end
        
    end
    if noShift==0, %last unsuccessful iteration will not be used...
        shiftPos = [shiftPos, lastShiftPos];
        shiftDirs = [shiftDirs, lastShiftDir];
    end
    if lastShiftDir == 0,
        shiftPos = []; shiftDirs = [];
    else
        meanBL = meanBL2;
    end
end

if ~isempty(shiftPos),
    for i = 1:length(shiftPos),
        if shiftDirs(i)==1,
            temp = find(bl_diag(max(1,shiftPos(i)-maxSearchLength):shiftPos(i))<bl_diag(shiftPos(i))*0.95,1,'last')+max(1,shiftPos(i)-maxSearchLength)-1;
            if isempty(temp)
                temp2 = find(bl_diag(max(1,shiftPos(i)-maxSearchLength2):shiftPos(i))<bl_diag(shiftPos(i))*0.975,1,'last')+max(1,shiftPos(i)-maxSearchLength2)-1;
                if isempty(temp2),
                    shiftStartPoint(i) = max(1,shiftPos(i) - maxSearchLength2);
                else
                    shiftStartPoint(i) = temp2;
                end
            else
                shiftStartPoint(i) = temp;
            end
            
            temp2 = find(bl_diag(shiftPos(i):min(shiftPos(i)+maxSearchLength,length(bl_diag)))>bl_diag(shiftPos(i))*1.05,1,'first')+shiftPos(i)-1;
            if isempty(temp2)
                temp3 = find(bl_diag(shiftPos(i):min(shiftPos(i)+maxSearchLength2,length(bl_diag)))>bl_diag(shiftPos(i))*1.025,1,'first')+shiftPos(i)-1;
                if isempty(temp3),
                    shiftEndPoint(i) = min(length(bl_diag),shiftPos(i)+maxSearchLength2);
                else
                    shiftEndPoint(i) = temp3;
                end
            else
                shiftEndPoint(i) = temp2;
            end
        else
            temp = find(bl_diag(max(1,shiftPos(i)-maxSearchLength):shiftPos(i))>bl_diag(shiftPos(i))*1.05,1,'last')+max(1,shiftPos(i)-maxSearchLength)-1;
            if isempty(temp)
                temp2 = find(bl_diag(max(1,shiftPos(i)-maxSearchLength2):shiftPos(i))>bl_diag(shiftPos(i))*1.025,1,'last')+max(1,shiftPos(i)-maxSearchLength2)-1;
                if isempty(temp2),
                    shiftStartPoint(i) = max(1,shiftPos(i) - maxSearchLength2);
                else
                    shiftStartPoint(i) = temp2;
                end
            else
                shiftStartPoint(i) = temp;
            end
            temp2 = find(bl_diag(shiftPos(i):min(shiftPos(i)+maxSearchLength,length(bl_diag)))<bl_diag(shiftPos(i))*0.95,1,'first')+shiftPos(i)-1;
            if isempty(temp2)
                temp3 = find(bl_diag(shiftPos(i):min(shiftPos(i)+maxSearchLength2,length(bl_diag)))<bl_diag(shiftPos(i))*0.975,1,'first')+shiftPos(i)-1;
                if isempty(temp3),
                    shiftEndPoint(i) = min(length(bl_diag),shiftPos(i)+maxSearchLength2);
                else
                    shiftEndPoint(i) = temp3;
                end
            else
                shiftEndPoint(i) = temp2;
            end
        end
    end
else
    shiftStartPoint = [];
    shiftEndPoint = [];
end

shiftAreaSig = ones(1,length(fhr));
for i = 1:length(shiftPos),
    shiftAreaSig(shiftStartPoint(i):shiftEndPoint(i)) = 0;
end;
shiftAreaSig2 = zeros(1,length(fhr));
for i = 1:length(shiftPos),
    shiftAreaSig2(shiftStartPoint(i):shiftEndPoint(i)) = 2;
end;
shifSigDiff = shiftAreaSig+shiftAreaSig2;

% Adjustment of the newBl_accDec to the newly set confidence intervals
for i=1:length(bl_accDec)
    if bl_accDec(i)> bl_diag(i)+shifSigDiff(i)*unitShiftDiff;
        bl_accDec(i) = bl_diag(i)+shifSigDiff(i)*unitShiftDiff;
    elseif bl_accDec(i)< bl_diag(i)-shifSigDiff(i)*unitShiftDiff;
        bl_accDec(i) =  bl_diag(i)-shifSigDiff(i)*unitShiftDiff;
    end
end

%%
function [shiftUpOut,lastShiftDir, meanBL2] = checkUp(bl_diag, meanBL, changeUp, segmentLength, fs)
shiftUpOut = []; lastShiftDir = +1;
meanBL2 = mean(bl_diag(changeUp:min(length(bl_diag),changeUp+segmentLength*60*fs)));
if  1.1*meanBL < meanBL2,
    shiftUpOut =  [shiftUpOut; changeUp];
end

%%
function [shiftDownOut,lastShiftDir, meanBL2] = checkDown(bl_diag, meanBL, changeDown, segmentLength, fs)
shiftDownOut = []; lastShiftDir = -1;
meanBL2 = mean(bl_diag(changeDown:min(length(bl_diag),changeDown+segmentLength*60*fs)));
if  0.9*meanBL > meanBL2,
    shiftDownOut =  [shiftDownOut; changeDown];
end

function [bradyVect, bradyPosOut] = eFIGO_detectBrady(fhr,fs,decel)
%%
bradyPosOut = [];
bradyAmpLimit = 110;
bradyTimeLimit = 2*60;%s

%%
counter = 0;
for i = 1:length(decel),
    if mean(fhr(decel(i).start:decel(i).end)) < bradyAmpLimit
        if (decel(i).end - decel(i).start)/fs > bradyTimeLimit
            counter = counter + 1;
            bradyPosOut(counter).start = decel(i).start;
            bradyPosOut(counter).end = decel(i).end;
            bradyPosOut(counter).iDecel = i;
        end
    end
end

%%
bradyVect = NaN(1,length(fhr));
if ~isempty(bradyPosOut),
    for i = 1:length(bradyPosOut),
        bradyVect(bradyPosOut(i).start:bradyPosOut(i).end)=1;
    end
end
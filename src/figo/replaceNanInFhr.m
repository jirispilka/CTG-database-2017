function nxData = replaceNanInFhr(nxData) 

%%% START OF TEMPORARY CODE TO MAKE IT WORK WITH CTU-UHB Database
addedZerosAtEnd = [];
if length(find(isnan(nxData)==1))>0 || length(find(nxData==-1))>0,
    if ~isempty(find(nxData==-1)),
        addedZerosAtEnd = length(find(nxData==-1));
        nxData(find(nxData==-1))=[];
    end
    lengthOfData = length(nxData);
    timeOldData = 1:1:lengthOfData;
    oldData = nxData;
    timeNewData = timeOldData;
    timeNewData(find(isnan(nxData)==1))=[];
    nxData(find(isnan(nxData)==1))=[];
    data_interp = pchip(timeNewData,nxData,timeOldData);
    nxData = data_interp;
    %figure; plot(data_interp,'r'); hold on; plot(oldData);
end
if ~isempty(addedZerosAtEnd)
    clear('nxData');
    nxData = [data_interp,ones(1,addedZerosAtEnd)];
end
if size(nxData,1)~=1,
    nxData = nxData';
end
%%% END OF TEMPORARY CODE TO MAKE IT WORK WITH CTU-UHB Database

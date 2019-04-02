%% Compute the features for all signal
clear;
close all; clc;

%% parameters need to be set

win_min = 20; % window length for analysis
winStep_min = 5; % time step

% interpolate small gas with length <= XX
nLimSmallGaps = 50;

% features to compute
sFeatSelect = 'stat';
% sFeatSelect = 'spectral';
% sFeatSelect = 'figo'; % for segments longer than 20 min.
% sFeatSelect = 'mfhurst'; % requires toolbox from Herwig, Patrice, Roberto

% example of combination
% sFeatSelect = 'figo_spectral_mfhurst';
% sFeatSelect = 'stat_figo_spectral_mfhurst';

bplot = 0;

sname_save = 'Features_CTU_stat_20190329.mat';
% sname_save = 'Features_CTU_stat_spectral_figo_mfhurst_20190329.mat';

%% data 
path_db = '~/data/CTU_UHB_2017/';
dirs = {'matfiles'};

aFiles = {};
for i = 1:length(dirs)   
    files = getAllFiles(fullfile(path_db,dirs{i}), 'mat');
    aFiles = [aFiles; files];
end

%% parameters (remains fixed)

fs = 4;
nNr = length(aFiles);

nWin_samp = win_min*60*fs;
nWin_limit = nWin_samp/2;

nWinStep_samp = winStep_min*60*fs;
nOverlap_samp = nWin_samp - nWinStep_samp;

%% init
pH = zeros(nNr,1);
acFeatures = cell(nNr, 1);
acSegmentInfo = cell(nNr, 1);
acMetainfo = cell(nNr, 1);

%% process files
cnt = 0;

for i = 1:nNr
    %tic
    name = aFiles{i};
    fprintf('%d/%d, name: %s \n',i,nNr,name);
       
    c = load(name);
    %data = c.bpm;
    data = c.bpm_nan;
    
    data = interpolateSmallGaps(data,fs,bplot,nLimSmallGaps);
    %data = interpolateSmallGaps(data,fs,bPlot,10e10);
    
    %% info
    pH(i) = c.info.pH;
    
    ind_stageII = fix(c.info.ind_stageII);
    ind_stageII = min(ind_stageII, length(data));
    
    %% Align the data by position of stage II
    % prepend NaN to stage I to fit the sliding window
    nr_steps = ceil(length(data(1:ind_stageII))/nWinStep_samp);
    nr_prepend = abs(ind_stageII - nr_steps * nWinStep_samp) + 1;

    % prepend NaN
    if nr_prepend ~= 0
        nanpad = nan(nr_prepend,1);
        data = [nanpad; data];
        ind_stageII = nr_prepend + ind_stageII;
        ind_stageII = min(ind_stageII, length(data));
        beginI = 1;
        endI = ind_stageII - 1;
    end
    
    % if the second stage is not available
    if ind_stageII == length(data)
        dataI = data(beginI:endI);
        dataII = [];
    else
        % second stage is available
        dataI = data(beginI:endI);
        dataII = data(ind_stageII:end);
           
        %stageII_fhr_len_samp(i) = length(dataII);
        %stageII_fhr_len_min(i) = length(dataII)/fs/60;
        
        % if the second stage is shorter than XXX minutes
        nr_steps = ceil(length(dataII)/nWinStep_samp);
        r2 = round(nr_steps * nWinStep_samp - length(dataII));
        if r2 ~= 0
            nanpad = nan(r2,1);
            dataII = [dataII; nanpad];
        end
    end
    
    if bplot
        xTime = (1/fs:1/fs:length(data)/fs)/60;
        figure(10000+i); clf;
        hold on;
        plot(xTime,data,'k')
        plot(xTime(ind_stageII:end),data(ind_stageII:end),'r')
        plot(xTime(beginI:endI),data(beginI:endI),'b')
        grid on;
        legend('FHR','IIstage for analysis','Istage for analysis','Location','best')
    end
       
    data = [dataI;dataII];
    
    % Mark the segments with aSegStage
    % 1 = Stage I
    % 2 = Stage II
    % 12 = overalap between the stages
    r = nWin_samp - nOverlap_samp;
    kall = (length(data) - nOverlap_samp)/r; % number of segments
    kI = (length(dataI) - nOverlap_samp)/r;
    kII = (length(dataII) - nOverlap_samp)/r;
    
    aSegStage = 12*ones(1,kall);
    aSegStage(1:kI) = 1;
    ind = length(aSegStage) - kII + 1;
    aSegStage(ind:end) = 2;
    
    % segments are indexed from the end
    % aSegStageI_index - index segments from stageI
    % aSegStageII_index - index segments from stageII
    aSeg_index = length(aSegStage):-1:1;
    aSegStageI_index = zeros(1,kall);   
    aSegStageI_index(1:kI) = kI:-1:1;
    
    aSegStageII_index = zeros(1,kall);
    if kII > 0
        aSegStageII_index(ind:end) = kII:-1:1;
    end
    
    %% compute on whole data (FIGO)
    if contains(sFeatSelect,'figo') || strcmpi(sFeatSelect,'all')             
        % remove NaN
        if i == 1709 
            data(33239) = NaN;
        end
        [temp,aGapBeg,aGapEnd] = removeNaNsAtBeginAndEnd(data);
        temp = [nan(size(aGapBeg))'; temp];
        
        [eFigo] = analyzeFHR_enhancedFIGO(temp', fs, '', bplot);
        % add NaN back
        aGapEnd = nan(size(aGapEnd));
        eFigo.baseLine = [eFigo.baseLine,aGapEnd];
        eFigo.baseLineAccDecc = [eFigo.baseLineAccDecc,aGapEnd];
        eFigo.multVect = [eFigo.multVect,aGapEnd];
        eFigo.decels.vector = [eFigo.decels.vector,aGapEnd];
        eFigo.accels.vector = [eFigo.accels.vector,aGapEnd];
        eFigo.brady.vector = [eFigo.brady.vector,aGapEnd];
        % compute the features
        featfigo = uti_computeFigoSegment(data,[1,length(data)],eFigo,fs);
        featfigo.bslnAllBeta0 = featfigo.bslnBeta0;
        featfigo.bslnAllBeta1 = featfigo.bslnBeta1;
        eFigo = copyStruct(featfigo,eFigo);

    else
        eFigo = [];
    end
    
    %% prepare segments
    
    [aSeg,aaBegEnd,~] = extractSegments(data, nWin_samp, nOverlap_samp);
    cTempFeat = cell(size(aSeg,2),1);

    % feat - NaN prototype
    if ~exist('featNaN','var')
        featNaN = uti_computeFeatPreSelected(nan(nWin_samp,1),sFeatSelect,aaBegEnd(1,:),fs,eFigo);
    end
    
    %% Compute features on segments
    aquality = zeros(size(aSeg,2),1);
    for k = 1:size(aSeg,2)
        %fprintf(' segment %d/%d\n',k,size(aSeg,2));
        seg = aSeg(:,k);
        bNaN = false;
                
        aquality(k) = 1-sum(isnan(seg))/length(seg);
        if sum(isnan(seg)) > nWin_limit
            bNaN = true;
        end
        
        if aSegStage(k) > 1
            [~, ~, gap_end] = removeNaNsAtBeginAndEnd(seg);
            if length(gap_end) >= nWinStep_samp
                bNaN = true;
            end
        end
        
        % feature computation if not NaN
        if bNaN
            cTempFeat{k} = featNaN;
        else
            cTempFeat{k} = uti_computeFeatPreSelected(seg,sFeatSelect,aaBegEnd(k,:),fs,eFigo);
        end
    end
    
    % fill in segment information
    cTempSeg = cell(size(aSeg,2),1);
    for k = 1:size(aSeg,2)
        val = strsplit(c.info.fileNameMat, '_');
        c_seg.name = c.info.dbID;
        c_seg.pH = c.info.pH;
        c_seg.year = int16(str2double(val(2)));
        c_seg.segStart_samp = aaBegEnd(k,1);
        c_seg.segEnd_samp = aaBegEnd(k,2);
        c_seg.segStage = aSegStage(k);
        c_seg.segIndex = aSeg_index(k);
        c_seg.segStageI_index = aSegStageI_index(k);
        c_seg.segStageII_index = aSegStageII_index(k);
        cTempSeg{k} = c_seg;
        clear c_seg
    end
    
    % reformat matrix
    [aFeatMatx, aFeatNames] = createFeatMatxFromCell(cTempFeat);
    [aSegMatx, aSegNames] = createFeatMatxFromCell(cTempSeg);

    acFeatures{i} = aFeatMatx;
    acSegmentInfo{i} = aSegMatx;
    acMetainfo{i} = c.info;

end

%% save to matlab file
if exist('sname_save','var')
    save(sname_save, 'aFeatNames', 'aFiles', 'sFeatSelect', ...
        'nWin_samp','nWin_limit', 'winStep_min','win_min', ...
        'acFeatures', 'acSegmentInfo', 'acMetainfo');
end

%% create csv file

sname_save_csv = strrep(sname_save, 'mat', 'csv');

% % write header
header = [aSegNames, aFeatNames];
fid = fopen(sname_save_csv, 'w+');
for i = 1:length(header)
    fprintf(fid, '%s', header{i});
    if i ~= length(header)
        fprintf(fid, ',');
    end
end
fprintf(fid, '\n');

aNames = [aSegNames, aFeatNames];
for i = 1:length(acFeatures)
    cfeat = acFeatures{i};
    cseg = acSegmentInfo{i};
    dlmwrite(sname_save_csv, [cseg, cfeat],'-append');
    
    t = array2table([cseg, cfeat], 'VariableNames', aNames);
    if i == 1
        T = t;
    else
        T = [T; t];
    end
end

T.name = int32(T.name);
writetable(T, sname_save_csv);



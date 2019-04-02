function cFeatAll = copyStruct(cFeat,cFeatAll)
% copy one struct to another
%
% Jiri Spilka
% Czech Technical University in Prague, 2015

if ~exist('cFeatAll','var')
    cFeatAll = [];
end

sNames = fieldnames(cFeat);
for i = 1:length(sNames)
    cFeatAll.(sNames{i}) = cFeat.(sNames{i});
end
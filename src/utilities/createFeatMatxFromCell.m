function [aFeatMatx, aFeaturesNames] = createFeatMatxFromCell(acFeatures)
% CREATEFEATMATXFROmCELL transform cell to an array
%
% Synopsis:
%   [aFeatMatx aFeaturesNames] = createFeatMatxFromCell(acFeatures)
%
% Description:
%   Transform cell to an array.
%   
% Input:
%  cFeatures [nx1 cell]
% 
% Output:
%  aFeatMatx [nxm]
%  aFeaturesNames [nxstring]
%
% Examples:
%  createFeatMatxFromCell(acFeatures)
%
% See also: 
%
% About: 
% Jiri Spilka
% Czech Technical University in Prague, 2015
%
% Modifications:
%

nNrRecords = size(acFeatures,1);

% because some of the cell could be empty
for i = 1:length(acFeatures)
    if ~isempty(acFeatures{i})
        idxNotEmpty = i;
        break;
    end
end

aFeaturesNames = fieldnames(acFeatures{idxNotEmpty})';
nNrFeatures = length(aFeaturesNames);

aFeatMatx = zeros(nNrRecords, nNrFeatures);

for i = 1:nNrRecords
    cFeatTemp = acFeatures{i};
    
    if isempty(cFeatTemp)
        continue;
    end
    
    for k = 1:nNrFeatures
        sName = aFeaturesNames{k};
        aFeatMatx(i,k) = cFeatTemp.(sName);
    end
end
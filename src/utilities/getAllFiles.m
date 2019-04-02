function fileList = getAllFiles(dirName,sExtension,sMask)
% GETALLFILES get all files within directory
%
% Synopsis:
%  [fileList] = getAllFiles(fname,[sExtension,[sMask]])
%
% Description:
%  This function finds all files within a directory. A directory is
%  recursively searched to subdirectories. Founed files are filtered by on
%  extension or by a mask
%
% Input:
%  dirName - directory name
%  sExtension - [optional] file extension, including dot
%  sMask - [optional] filter mask
% 
% Output:
%  fileList - [nxcell] founded files satisfying the conditions
%
% Example:
%   fileList = getAllFiles('c:\data'); 
%   fileList = getAllFiles('c:\data','txt'); 
%   fileList = getAllFiles('c:\data','gz','bdi'); 
%
% See also 
%
% About: 
%  Jiri Spilka
%  Czech Technical University in Prague
%  http://people.ciirc.cvut.cz/~spilkjir
%
% Modifications:
% 

if nargin == 1
    sExtension = '';
    sMask = '';
elseif nargin == 2
    sMask = '';
end

dirData = dir(dirName);
% Find the index for directories
dirIndex = [dirData.isdir];
sMask = lower(sMask);

% Get a list of the files
fileList = {dirData(~dirIndex).name}';
% filter files by an extension
fileList = getAllFilesFilterMask(fileList, sExtension,sMask);

% Prepend path to files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),fileList,'UniformOutput',false);
end

% Get a list of the subdirectories
subDirs = {dirData(dirIndex).name};

% Find index of subdirectories that are not '.' or '..'
validIndex = ~ismember(subDirs,{'.','..','.svn'});
for iDir = find(validIndex)
    nextDir = fullfile(dirName,subDirs{iDir});
    % Recursively call getAllFiles
    fileList = [fileList; getAllFiles(nextDir,sExtension,sMask)];
end

%% applyFilterMask - returns files with the extension
function filesOut = getAllFilesFilterMask(files, sExtension,sMask)

N = length(files);
aIndciesAccept = zeros(N,1);
for i = 1:N
    if (fileEndsWith(files{i},sExtension))
        if isempty(sMask) || ~isempty(strfind(lower(files{i}),sMask))
            aIndciesAccept(i) = i;
        end
    end
end

filesOut = files(aIndciesAccept ~= 0);
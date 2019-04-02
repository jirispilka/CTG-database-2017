function b = fileEndsWith(file,sInputExt)
% FILEENDSWITH check if a file ends with predifined extension
%
% Synopsis:
%  [cSignal] = FileEndsWith(file,sInputExt)
%
% Description:
%  Returns true if a file ends with an extension
%
% Input:
%  file - [str] file name
%  sInputExt - [str] extension
% 
% Output:
%  b - [bolean] 
%
% Example:
%  b = FileEndsWith('text.txt','txt')
%
% See also: 
%
% About: 
% Jiri Spilka
% Czech Technical University in Prague, 2015
%
% Modifications:
%

if ~isempty(sInputExt)
    if (~strcmp(sInputExt(1),'.'))
        %sInputExt = strcat('.',sInputExt);
        sInputExt = ['.' sInputExt];
    end
    
    % Look for EXTENSION part
    ind = find(file == '.', 1, 'last');
    
    if isempty(ind)
        b = false;
        return;
    else
        ext = file(ind:end);
    end
    
    b = strcmp(ext,sInputExt);
else   
    b = true;
end
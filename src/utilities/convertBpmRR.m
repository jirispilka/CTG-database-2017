function [aOutData, idxAtf] = convertBpmRR(aInData)
% CONVERTBPMRR - convert the both ways: bpm2rr and rr2Bpm
%
% Jiri Spilka
% Czech Technical University in Prague, 2015

aOutData = 60*1000./aInData;
idxAtf = isinf(aOutData);
aOutData(idxAtf) = 0;
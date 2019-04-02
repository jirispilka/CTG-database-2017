function [SSI,result_out,FHR,bandwidth,density,xmesh] = baselineKDE_comp(nxData,aprioriConf,n,nTrshUnstable,nTrshModStable)
% Method using kernel density estimate (KDE) for detection of baseline.
% 
% INPUTS:
%       nxData  - vector of data from which the density estimate is constructed;
%       n     - the number of mesh points used in the uniform discretization; 
%               n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12, default 256;
%       unst  - treshold for unstable (SSI <= unst), default 0.02;
%       mst   - treshold for moderately stable (unst < SSI <= mst), default
%               0.05;
%
% OUTPUTS:
%       SSI         - signal stability index, maximum of density;
%       result      - function's result, treshold: <= 0.02 unstable equals "3",
%                     <= 0.05 moderately stable equals "2", > 0.05 stable equals "1";
%       FHR         - baseline value for SSI point;
%       bandwidth   - the optimal bandwidth (Gaussian kernel assumed);
%       density     - column vector of length 'n' with the values of the density
%                     estimate at the grid points;
%       xmesh       - the grid over which the density estimate is computed;
%
%      If no output is requested, then the code automatically plots a graph
%      of the density estimate.
%
% Authors:
%       Lukas Zach, Vaclav Chudacek 
%
% Updated:
%       Vaclav Chudacek 30/11/2016
%
% Reference: 
%       Title:      Computerized fetal heart rate analysis in labor: 
%                   detection of intervals with un-assignable baseline
%       Authors:    Antoniya Georgieva et al.
%       LibName:    Georgieva2011
%% Input arguments revised
if nargin < 3
    n = 64;
elseif nargin < 4
    nTrshUnstable = 0.02; 
    nTrshModStable = 0.05; 
elseif nargin < 5
    nTrshModStable = 0.05; 
end

%% Main call
[bandwidth,density,xmesh] = kde(nxData,n);
[SSI,k] = max(density);
FHR = xmesh(k);

if SSI <= nTrshUnstable || FHR < 50 || aprioriConf == -1,
	result_out = 3; 
elseif SSI <= nTrshModStable
	result_out = 2;
else
	result_out = 1;
end


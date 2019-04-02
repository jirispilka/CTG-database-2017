function [feat,esta,logstat] = featuresMF_RFL(data,j1,j2,gamint,pp)
% PA Lyon, May 2013
% gathers various MF Leaders versions
% Structure Function
% Quantile
% Bayesian

% code from Roberto - can handle the missing data

if ~exist('pp','var')
    pp = 0; % pp = 0; % norm of p-leaders
end

%% >---- NEW PARAMETERS ---->
% p-Leaders
% pp=1; % norm of p-leaders
      % if pp=0: compute standard leaders
% initial scaling range
J1LF=1;
% symmetrize wavelet (1d only)
sym=0;              
% Select use of scaling correction term:
NLS = 2;          % 0: don't use correction
                  % 1: use correction if $p \neq \infty$
                  % 2: use correction if $p \neq \infty$ and $\eta(p) \geq 0$

% check and write parameters:
%checkParam_MF_BS_tool;

%%
N = length(data);

%% Analysis Params

%--- Analysis Methods to use
methodWT=[1,2];   % 1 - DWT  (discrete wavelet transform coefficients)
                  % 2 - PWT  (p-Leaders)
                  
%% parameters for different exponents
dgam=+0.05;          % differential fractional integration (+) or derivation (-) for estimation oscillation and cancellation exponent
deltap=+0.02;       % delta 1/ for derivative d/d(1/p) 
Texp=0;             % 0: p-exponent.            
                    % 1: oscillation exponent.      [dgam used]
                    % 2: lacunary exponent.         [deltap used]
                    % 3: cancellation exponent.     [dgam used]
                    %       NOTE: need gamint=0
                    
%%
%--- Estimation Parameters
                  
MomNul=3;           % vanishing moments wavelet (Daubechies')

% gamint=0;         % Wavelet domain fractional integration by gamint

% j1=3;             % j1,j2 - scaling range
% j2=15;

wtype=0;            % linear regression : 0 - ordinary ; 1 - weighted nj(j) ; 2 - weighted with estimated variance
Fun=111;            %       number xyz - what is to be calculated:  
                    %       x : zeta(q) [0 or 1]    (scaling exponents)
                    %       y : D(h) [0 or 1]       (multifractal spectrum)
                    %       z : cp [0 or 1]         (log-cumulants)
                    %       e.g. param.EstFun=001 calculates only Cp, 
                    %            param.EstFun=110 calculates SF zeta(q) and D(h)

Cum=4;              % Highest order of log cumulant cp to estimate

nq = 12 ; qm = 0.1 ; qM = 10 ;
q = logspace(log10(qm),log10(qM),nq) ;
q = [-fliplr(q), -2, -1 0, 1, 2 q] ;
q = sort(q) ;
nq = length(q) ;

%--- BOOTSTRAP
% Bootstrap parameters
T_S=0;              % time - scale block [1] or time block [0] bootstrap
B1=0;              % # primary resamples - NO BOOTSTRAP (ESTIMATION ONLY): B1=0;
B2=0;               % # double bootstrap resamples
CI=0;               % calculate confidence intervals
TEST=0;             % calculate tests
if T_S              % Block Length
    Block=floor(N/32);     
else
    Block=2*MomNul; 
end
Method=[3];         % Methods to be used [1-Normal; 2-Basic; 3-Percentile; 4-Studentised; 5-Adjusted Basic; 6-Adjusted Percentile]
Alpha=0.1;          % Significance Level (1-Alpha)
Jflag=1;            % 1: do not do Bootstrap on scales that are not implied in regression (use 0 only if you need confidence intervals for scales not involved in regressions)
% Test Parameters
Tnull=zeros(1,Cum); % Null Hypothesis Parameter (for cp only)
Type=[4];           % Null Hypothesis H0: T=Tnull Types [  HA: T>Tnull;   HA: T<Tnull;   H0: |T-Tnull|=0;   H0: T-Tnull=0  ]

%--- VERBOSITY PARAMETERS
FigNum=0;           % 0 - no Figures
verbose=0;          % Verbosity level:
                    %   0 : batch mode
                    %   1 : estimate, CI, Test table
                    %       Figure estimate if FigNum
                    %   11: Figure estimate if FigNum
                    %   2 : 1 + Figure LogScale if FigNum
                    %   21: 11+ Figure LogScale if FigNum
                    %   3 : 2 + interactive mode (change Estimates, j1, j2, wtype,
                    %           Analysis method, Alpha, Text output, Figures)

levelsets = [0.9 0.7 0.5 0.3] ;
nlevelsets = length(levelsets) ;

                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));
index = find(q==2);
% check and write parameters:
checkParam_MF_BS_tool_singul;
[paramEST, paramBS, paramTest]=MF_BS_tool_param_singul(MomNul,gamint,dgam,deltap,Texp,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S,CI,TEST,pp,J1LF,sym,NLS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
[esta,conf,signif,logstat] = MF_BS_tool_singul(data, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
% toc;

% %% OLD CODE
% % check and write parameters:
% checkParam_MF_BS_tool_NLpL_LF;
% [paramEST, paramBS, paramTest]=MF_BS_tool_param_NLpL_LF(MomNul,gamint,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S,CI,TEST,pp,J1LF,sym);
% [esta,conf,signif,logstat] = MF_BS_tool_NLpL_LF(data, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);

%% classical analysis with predefined scales
j1min = logstat.j1min; 
j2max = logstat.j2max;

C1_j = logstat.LWT.est(end-Cum+1,j1min:j2max);
C2_j = logstat.LWT.est(end-Cum+2,j1min:j2max);
C3_j = logstat.LWT.est(end-Cum+3,j1min:j2max);
C4_j = logstat.LWT.est(end-Cum+4,j1min:j2max);

C1 = esta.LWT.t(end-Cum+1) ;
C2 = esta.LWT.t(end-Cum+2) ;
C3 = esta.LWT.t(end-Cum+3) ;
C4 = esta.LWT.t(end-Cum+4) ;
hmin_noint = esta.DWT.h_minnoint ;
zetaqL(1:length(q)) = esta.LWT.t(1:length(q)) ;
zetaqD(1:length(q)) = esta.DWT.t(1:length(q)) ;
H_DWT = zetaqD(index)/2 - gamint;
H_LWT = zetaqL(index)/2 - gamint;

estid = find(logstat.q==2);
scales=log2(logstat.DWT.scale);
temp=logstat.DWT.est(estid,:);
H_j = temp(j1min:j2max);

% figure(1025); clf;
% hold on;
% grid on;
% plot(scales, temp, 'k.-');

%% SCALING RANGE DETECTION
% 
% jlimits = [3 6 ; 8 11];
% %-------------------------------------------------------------------------------
% % Detection for H from Wavelet Coeff:
% qwanted = 2;
% estid = find(paramEST.q == qwanted);
% 
% logsf = logstat.DWT.est(estid,:);   % Selected estimates, either DWT or LWT
% logsfBS = logstat.DWT.estB(:,estid,:);  % Selected boostrap resamples
% logsfBBS = [];    % Do not use doble bootstrap
% 
% [j1H_D,j2H_D,gof, slope, intercept, Q] = rangeSelQ_v14(logsf, logsfBS, logsfBBS, jlimits);
% if verbose
%     fprintf( 'Detected scaling range for H_DWT: ( %g, %g )\n', j1H_D,j2H_D)
% end
% H_DWTadt = slope/qwanted - gamint;
% 
% %-------------------------------------------------------------------------------
% % Detection for c1 only:
% cwanted = 1;   % Select c1
% estid = length(paramEST.q) * 3 + cwanted;
% logsf = logstat.LWT.est(estid,:);   % Selected estimates, either DWT or LWT
% logsfBS = logstat.LWT.estB(:,estid,:);  % Selected boostrap resamples
% logsfBBS = [];    % Do not use doble bootstrap
% 
% [j1c1,j2c1,gof, slope, intercept, Q] = rangeSelQ_v14(logsf, logsfBS, logsfBBS, jlimits);
% if verbose
%     fprintf( 'Detected scaling range for c1: ( %g, %g )\n', j1c1,j2c1 )
% end
% C1adt = log2(exp(1))*slope;
% esta.LWT.t(end-1:end)
% 
% %-------------------------------------------------------------------------------
% % Detection for c2 only:
% cwanted = 2;   % Select c2
% estid = length(paramEST.q) * 3 + cwanted;
% logsf = logstat.LWT.est(estid,:);   % Selected estimates, either DWT or LWT
% logsfBS = logstat.LWT.estB(:,estid,:);  % Selected boostrap resamples
% logsfBBS = [];    % Do not use doble bootstrap
% 
% [j1c2,j2c2,gof, slope, intercept, Q] = rangeSelQ_v14(logsf, logsfBS, logsfBBS, jlimits);
% if verbose
%     fprintf( 'Detected scaling range for c2: ( %g, %g )\n', j1c2,j2c2 )
% end
% C2adt = log2(exp(1))*slope;
% % % esta.LWT.t(end-1:end)
% % 
% % %-------------------------------------------------------------------------------
% % Detection for c1 and c2 jointly:
% cwanted = 1:2;   % Select c1 and c2 jointly
% estid = length(paramEST.q) * 3 + cwanted;
% logsf = logstat.LWT.est(estid,:);   % Selected estimates, either DWT or LWT
% logsfBS = logstat.LWT.estB(:,estid,:);  % Selected boostrap resamples
% logsfBBS = [];    % Do not use doble bootstrap
% 
% [j1c1c2,j2c1c2,gof, slope, intercept, Q] = rangeSelQ_v14(logsf, logsfBS, logsfBBS, jlimits);
% if verbose
%     fprintf( 'Detected scaling range for c1 and c2 jointly: ( %g, %g )\n', j1c1c2,j2c1c2 )
% end
% C1adt_c1c2 = log2(exp(1))*slope(1);
% C2adt_c1c2 = log2(exp(1))*slope(2);
% % esta.LWT.t(end-1:end)

%%
% % Leaders quantiles
% [coef, leaders, nj] = DxPx1d_LF(data, MomNul, gamint) ;
% %       j2 = min(j2,length(leaders)) ;
% nqt=25; qt=linspace(0,1,nqt+1);
% for j=1:j2;
%     if j > length(leaders)
%         continue;
%     end
%     lx=sort(log2(leaders(j).value), 'ascend');nlx=length(lx);
%     jqt=round((nlx/2).^qt);jqt=[jqt nlx-fliplr(jqt(2:end))]; jqt=max(1,min(nlx,jqt));
%     Qj(:,j)=lx(jqt);
% end
% yj=Qj; varyj=ones(size(yj)); 
% [hQt]=MFA_BS_regrmat(yj,varyj,nj.L,wtype, j1,j2);
% dqt=[qt fliplr(qt(1:end-1))];
% indexmedian = find(dqt==1) ;
% C1median = squeeze(hQt(indexmedian)) ;
% 
% for z=1:1:nlevelsets
%     level = levelsets(z) ;
%     indexlevel = find(dqt>=level) ;
%     index1 = indexlevel(1) ;
%     index2 = indexlevel(end) ;
%     Dp(z) = squeeze(hQt(index1))-C1median ;
%     Dm(z) = -(squeeze(hQt(index2))-C1median) ;
%     Dpm(z) = squeeze(hQt(index1)-hQt(index2)) ;
% end
% z = 1 ;
% C2median = -squeeze(Dpm(z)) ;
% C2medianp1 = squeeze(Dp(z)) ;
% C2medianm1 = squeeze(Dm(z)) ;
% C2medianpm1 = squeeze(Dp(z))./squeeze(Dm(z)) ;
% z = 2 ;
% C2median2 = -squeeze(Dpm(z)) ;
% C2medianp2 = squeeze(Dp(z)) ;
% C2medianm2 = squeeze(Dm(z)) ;
% C2medianpm2 = squeeze(Dp(z))./squeeze(Dm(z)) ;

%%
if gamint == 0
    feat.MF_HDWT = H_DWT;
    %feat.MF_HLWT = H_LWT;

    feat.MF_hmin_noint = hmin_noint;    
else
    
    feat.MF_HDWT = H_DWT;
    %feat.MF_HLWT = H_LWT; 
    feat.MF_hmin_noint = hmin_noint;
    
    feat.MF_c1 = C1;
    feat.MF_c2 = C2;
    feat.MF_c3 = C3;
    feat.MF_c4 = C4;
    
    %feat.c1_j = C1_j;
    %feat.c2_j = C2_j;
    %feat.c3_j = C3_j;
    %feat.c4_j = C4_j; 
    %feat.H_j = H_j;
    
%     feat.MF_HDWT_adt = H_DWTadt;
%     feat.MF_HDWTj1 = j1H_D;
%     feat.MF_HDWTj2 = j2H_D;       
%     
%     feat.MF_c1adt = C1adt;
%     feat.MF_c1j1 = j1c1;
%     feat.MF_c1j2 = j2c1;    
%     
%     feat.MF_c2adt = C2adt;
%     feat.MF_c2j1 = j1c2;
%     feat.MF_c2j2 = j2c2;
%     
%     feat.MF_c1adt_c1c2 = C1adt_c1c2;
%     feat.MF_c2adt_c1c2 = C2adt_c1c2;
%     feat.MF_c1c2j1 = j1c1c2;
%     feat.MF_c1c2j2 = j2c1c2; 
% 
%     feat.C1median = C1median;
%     feat.C2median = C2median;
end
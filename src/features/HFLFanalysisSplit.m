function [LFHFratio,LFpow,HFpow,totPow,PowSplit,alpha] = HFLFanalysisSplit(sig,Fs,LFband,HFband,h, nVerbose)
% Routine for doing HF/LF analysis
%
% Usage
%   [LFHFratio,LFpow,HFpow,totPow] = HFLFanalysis(sig,Fs,LFband,HFband,h)
% Inputs
%   sig     signal to analyze
%   Fs      sampling frequency in Hz
%   LFband  vector defining low-frequency band in Hz
%           Default: [0.04 0.15]
%   HFband  vector defining low-frequency band in Hz
%           Default: [0.15 ]
%   h       the spectral estimator (type `help spectrum' to see what
%           options are available
%           Default: Welch spectral estimator with segment length 1024
%
% Outputs
%   LFHFratio   the ratio `power in the LF band'/`power in the HF band'
%   LFpow       `power in the LF band'/`total power in the HF and LF bands'
%   HFpow       `power in the HF band'/`total power in the HF and LF bands'
%   totPow      total power in the HF and LF bands
%
% Copyright (c) 2010, Hannes Helgason, Patrice Abry, Jiri Spilka

% Modifications
% 2013-11-06: JS - used trapeizod integration instead of rectangular
% 2013-11-06: JS - interpolate the spectrum to better find LF, HF boundaries
% 2013-11-05: JS - use function pwelch instead of psd (not working on my matlab)

if nargin < 4
    HFband = [0.15 0.4]; % in Hz (the previous default was [0.15 1])  
    if nargin < 3
        LFband = [0.04 0.15]; % in Hz
    end
end

%nfft = length(sig);
%f = (Fs/2)/nfft*(0:nfft-1); % frequency vector

% Instantiate spectrum object and call its PSD method.
if nargin < 5 || isempty(h)
    h = spectrum.welch;
    h.SegmentLength = 1024;
end
%h = spectrum.periodogram; % uses rectangular window by default
%h.WindowName = 'kaiser';

if nargin < 6
    nVerbose = 0;
end

% subtract the mean
sig = sig - mean(sig);

% % JS commented out, not working on my Matlab 2012a
% sigPsd= psd(h,sig,'Fs',Fs);
% f = sigPsd.Frequencies;
% P = sigPsd.Data;

% JS: 2013-11-05 - compute PSD using welch method, the psd was not working (Matlab 2012a) 
win = h.SegmentLength;
overlap = round(h.SegmentLength*0.8);%h.SegmentLength/2;
[P,f] = pwelch(sig,win,overlap,[],Fs,'onesided');

% interpolate the spectrum to better find frequency bands
% just for the integration
f_i = linspace(f(1),f(end), 2^4*length(f));
P_i = interp1(f,P,f_i,'linear');

% % %%%%%%%%%% rectangular interpolation %%%%%%%
% % % change frequency band to discrete frequencies
% f0 = LFband(1);
% f1 = LFband(end);
% ind = (f<=f1 & f>=f0);
% LFpow = sum(P(ind));
% 
% f0 = HFband(1);
% f1 = HFband(end);
% ind = (f<=f1 & f>=f0);
% HFpow = sum(P(ind));
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% trapeizodal interpolation %%%%%%%
LFpow = computeInbandEnergy(f_i,P_i,LFband(1),LFband(end));
HFpow = computeInbandEnergy(f_i,P_i,HFband(1),HFband(end));

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LFHFratio = LFpow/HFpow;
% normalize
totPow = HFpow+LFpow;
%HFpow = HFpow/totPow;
%LFpow = LFpow/totPow;

% spectral slope
ind = (0.04 <= f_i & f_i <= 1);
p2 = polyfit(log2(f_i(ind)),log2(P_i(ind)),1);
bestFit = polyval(p2,log2(f_i));

alpha = -p2(:,1);

%% plots

if nVerbose > 0
    fontsize = 14; ms = 4; lw = 2;

    figure(10258) ; clf 
    hold on;
    set(gca,'FontSize',fontsize) ; 
    
    plot(log2(f),log2(P), 'k','MarkerSize',ms,'LineWidth',lw); 
    plot(log2(f_i),log2(P_i), '--k','MarkerSize',ms,'LineWidth',1); 
    %stairs(log2(f),log2(P), 'r','MarkerSize',ms,'LineWidth',lw); 
    
    plot(log2(f_i),bestFit, '--r','LineWidth',1)    
    plot(log2(f_i(ind)),bestFit(ind), '-r','LineWidth',1)    
    
    a = axis;
    plot([log2(LFband(1)) log2(LFband(1))],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(LFband(end)) log2(LFband(end))],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(HFband(end)) log2(HFband(end))],[a(3) a(4)],'--b','LineWidth',2)
    xlabel('log(f)')
    ylabel('log(PSD)')
    legend('PSD')
    s =  sprintf('LFpow=%2.3f,HFpow=%2.3f,LF/HF=%2.2f',LFpow,HFpow,LFHFratio);
    title(s);
    grid on;    
end

%% Additional Bands PA Oct 2013 from TAsk Force paper
% European Heart Journal (1996) 17, 354-381

% rectangular interpolation
% ind = (f<0.003 & f>=0.0001);
% PowSplit(1) = sum(P(ind)) ;
% ind = (f<0.04 & f>=0.003);
% PowSplit(2) = sum(P(ind)) ;
% ind = (f<0.15 & f>=0.04);
% PowSplit(3) = sum(P(ind)) ;
% ind = (f<0.40 & f>=0.15);
% PowSplit(4) = sum(P(ind)) ;
% ind = (f<1.00 & f>=0.40);
% PowSplit(5) = sum(P(ind)) ;
% PowSplit(6) = sum(P(:)) ;

% trapeizodal interpolation
% the eps here is used to say lower than 0.003
PowSplit(1) = computeInbandEnergy(f_i,P_i,0.0001,0.003-eps); 
PowSplit(2) = computeInbandEnergy(f_i,P_i,0.003,0.04-eps);
PowSplit(3) = computeInbandEnergy(f_i,P_i,0.04,0.15-eps);
PowSplit(4) = computeInbandEnergy(f_i,P_i,0.15,0.40-eps);
PowSplit(5) = computeInbandEnergy(f_i,P_i,0.40,1.00-eps);
PowSplit(6) = computeInbandEnergy(f_i,P_i,0,0.04);

if nVerbose > 0
    h = figure(10259) ; clf 
    hold on;
    plot(log2(f_i),log2(P_i), 'k','MarkerSize',ms,'LineWidth',lw); 
    a = axis;
    plot([log2(0.0001) log2(0.0001)],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(0.003) log2(0.003)],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(0.04) log2(0.04)],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(0.15) log2(0.15)],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(0.40) log2(0.40)],[a(3) a(4)],'--b','LineWidth',2)
    plot([log2(1) log2(1)],[a(3) a(4)],'--b','LineWidth',2)
    text(-12.5,-8,'ULF','color','b') ;
    text(-8,-8,'VLF','color','b') ;
    text(-4,-8,'LF','color','b') ;
    text(-2.5,-8,'HF','color','b') ;
    text(-1,-8,'VHF','color','b') ;

    xlabel('log(f)')
    ylabel('log(PSD)')
    legend('PSD')
    s =  sprintf('Bands according to Task-Force\n');
    s =  sprintf('%s ULF: %2.1f, VLF: %2.1f, LF: %2.1f, HF: %2.1f, HHF: %2.1f \n tot= %2.1f', ...
        s,PowSplit(1),PowSplit(2),PowSplit(3),PowSplit(4),PowSplit(5),PowSplit(6));
    
    title(s);
    grid on;  
    fancyGraph(h);
end

%% compute band energy
function pow = computeInbandEnergy(f,Pxx,fl,fh)

ind = (fl <= f & f <= fh);
pow = trapz(f(ind),Pxx(ind));
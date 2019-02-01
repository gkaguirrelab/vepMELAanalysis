% function []=runVEPanalysis()
% A function that accesses other functions to run all the analysis code
% for Metropsis/VEP stimuli for individual observers.
%
% Syntax:
%  []=runVEPanalysis=()
%
% Description:
%	This function accesses other functions to run all the analysis code for Metropsis/VEP stimuli.
%
%
% Output: saves structure VEP
%   VEP                   - A structure that contains VEP data, TTL pulse,
%                           and timebase for both analog signals
%   audioRec              - A structure that contains mic data, and 
%                           sampling rate (Fs)
%   expParam              - contains observer ID, experiment ID, session 
%                           ID, data that was recorded from the VEP computer
%   mtrp                  - structure that contains all the metropsis data. 
%                           This includes observer ID, group (MWVA or HA free), 
%                           session, and temporal frequency for the stimuli 
%                           presented (TFtrials)
%   VDS                   - Visual discomfort scale values for the 36 trials

%% Load compiled data for a single observer
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');

filenameMAT=fullfile(getpref('vepMELAanalysis','melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID],['Exp' expID '_' observerID 'compiled.mat']);

open(filenameMAT);
VEP_main=ans.VEP;

%% Parse VEP data
[parsedVEPdata]=parseVEP(VEP_main);

%% Calculate TTF
vep_Fr=parsedVEPdata.vep_Fr;
TF=unique(parsedVEPdata.vep_Fr);
[ttf]=calcVEPttf(vep_Fr,parsedVEPdata.Fs,'TF',TF);

%% Plot TTF and VDS
% Plot mean Visual discomfort data

subplot(2,1,2)
VDSm=nanmean(parsedVEPdata.vds_Fr,2);
VDSstd=nanstd(parsedVEPdata.vds_Fr,[],2);
errorbar(A,VDSm,VDSstd,'-ok')
ylabel('visual discomfort scale')
xlabel('temporal frequency of stimulus')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 65];
ax.YLim=[0 10];

% Plot TTF
figure(5)
hold on
p_dataFr=squeeze(ttf.ttf);
P_Fm=mean(ttfFr,2);
P_Fstd=std(ttfFr,[],2);
errorbar(A,P_Fm,P_Fstd,'-ob')
title(expID)
ylabel('power spectra for stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 65];
ax.YLim=[0 0.02];
%end

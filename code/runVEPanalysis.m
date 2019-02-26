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
%                     c      session, and temporal frequency for the stimuli 
%                           presented (TFtrials)
%   VDS                   - Visual discomfort scale values for the 35 trials

%% Identify subject and create a path and file name to save data
observerID=input('observer ID:','s');
savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');
if ~exist(savePath,'dir')
    mkdir(savePath);
end

filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);
    
%% run all analyses for the 3 channel conditions
dur_in_sec=1.8;
starttime=0.2;

for x=1:3
    %% Load compiled data for a single observer from a single channel
    switch x
            case 1
                expID=input('experiment ID for luminance:','s');
                color='k';
            case 2
                expID=input('experiment ID for red/green:','s');
                color='r';
            case 3
                expID=input('experiment ID for blue/yellow:','s');
                color='b';
    end
    
    filenameMAT=fullfile(getpref('vepMELAanalysis','melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID],['Exp' expID '_' observerID 'compiled.mat']);

    open(filenameMAT);
    VEP_main=ans.VEP;

    %% Parse VEP data 
    [parsedVEPdata(x)]=parseVEP(VEP_main,'dur_in_sec',dur_in_sec,'starttime',starttime,'bandstop60',true,'bandstop120',true,'plot_all',true);
    Fs=VEP_main(1).vepDataStruct.params.frequencyInHz;
    XX=(1:length(parsedVEPdata(x).vep_Fr))/Fs;
    
    %% Calculate TTF
    [ttf(x)]=calcVEPttf(parsedVEPdata(x).vep_Fr,'normalize',true,'dur_in_freq',dur_in_sec*Fs,'plot_all',true);

    
    %% Plotting
    A=[1.625 3.25 7.5 15 30];
    % Plot mean Visual discomfort data
    figure(5)
    subplot(2,1,2)
    VDSm=nanmean(parsedVEPdata(x).vds_Fr,2);
    VDSstd=nanstd(parsedVEPdata(x).vds_Fr,[],2);
    errorbar(A,VDSm,VDSstd,['-o' color])
    hold on
    ylabel('visual discomfort scale')
    xlabel('temporal frequency of stimulus')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 65];
    ax.YLim=[0 10];

    
     % Plot TFF (power across averaged trials)
    
    figure(5)
    subplot(2,1,1)
    hold on
    errorbar(A,ttf(x).ttf_FrM,ttf(x).ttf_FrM-ttf(x).ttf_CI(:,1),ttf(x).ttf_CI(:,2)-ttf(x).ttf_FrM,['-o' color])
    title(observerID)
    ylabel('power spectra for stimulus frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 65];
    ax.YLim=[0 0.015];
    
    % Plot superimposed luminance, red/green, and blue/yellow in time
    % domain
    figure(11)
    for z=1:length(A)
        subplot(3,2,z)
        plot(XX,squeeze(mean(ttf(:,x).vep_Fr(z,:,:),2)),['-' color])
        title(['frequency=' num2str(A(z))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 parsedVEPdata(x).dur_in_freq/Fs];
        hold on
    end
    
    % Plot sum of power spectra averaged across repeat, and frequency for
    % each channel
    TTF_all=sum(mean(mean(ttf(x).ttf,2),1));
    
    figure(12)
    plot(x,TTF_all,['o' color])
    title('Sum TTF across frequencies');
    xlabel('Channel')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    hold on
    
%     % info to save per channel
%     VEP_main.mtrp.group;
%     VEP_main.mtrp.observer;

end



% save(filenameComp)
%end

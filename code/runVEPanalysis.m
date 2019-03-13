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
% Output: ***

%% Identify subject and create a path and file name to save data
observerID=input('observer ID:','s');
savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');
if ~exist(savePath,'dir')
    mkdir(savePath);
end

filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);
    
%% run all analyses for the 3 channel conditions
dur_in_sec=1.5;
starttime=0.5;

for x=1:3
    %% Load compiled data for a single observer from a single channel
    switch x
            case 1
                expID='LMS05';
                color='k';
                Color=[0.8 0.8 0.8];
            case 2
                expID='LM004';
                color='r';
                Color=[1 0.8 0.8];
            case 3
                expID='S0004';
                color='b';
                Color=[0.8 0.8 1];
    end
    
    filenameMAT=fullfile(getpref('vepMELAanalysis','melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID],['Exp' expID '_' observerID 'compiled.mat']);

    open(filenameMAT);
    VEP_main=ans.VEP;

    %% Parse VEP data
    
    [parsedVEPdata(x)]=parseVEP(VEP_main,'dur_in_sec',dur_in_sec,'starttime',starttime,'bandstop60',true,'bandstop120',true,'plot_sessions',false);
    Fs=VEP_main(1).vepDataStruct.params.frequencyInHz;
    XX=(1:length(parsedVEPdata(x).vep_Fr))/Fs;
    A=unique(VEP_main(x).mtrp.TFtrials);
    
    %% Process VEP data (gets rid of poor quality trials, and normalizes signal)
    [processedVEPdata(x)]=preprocessVEP(parsedVEPdata(x).vep_Fr, parsedVEPdata(x).vep_bkgd,'dur_in_sec',dur_in_sec,'normalize1',false);
    
    %% Calculate TTF
    [ttf(x)]=calcVEPttf(processedVEPdata(x).vep_Fr,'dur_in_sec',dur_in_sec,'plot_all',false,'TemporalFrequency',A);

    
    %% Plotting
    % Plot mean Visual discomfort data
    figure(5)
    subplot(2,1,2)
    VDSm=nanmedian(parsedVEPdata(x).vds_Fr,2);
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
    
    % Plot TTF
    if x==1
        figTTF=figure('Name',observerID);
    end
    figure(figTTF)
    hold on
    for YY=1:length(A)
        YYY=x+((YY-1)*3);
        subplot(5,3,YYY)
        TEMP=fill(cat(2,ttf(x).f(2:end),fliplr(ttf(x).f(2:end))),cat(2,squeeze(ttf(x).ttf_CI(YY,2:end,1)),fliplr(squeeze(ttf(x).ttf_CI(YY,2:end,2)))),Color,'EdgeColor','none');
        hold on
        plot(ttf(x).f,ttf(x).ttf_M(YY,:),['-' color])
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.YScale='log';
        ax.XLim=[min(ttf(x).f) 100];
        ax.XTick=A;
        ax.YLim=[0.00001 0.015];
        title(num2str(A(YY)))
        if YYY==7
            ylabel('Power')
        end
        
        if YYY==14
            xlabel('Frequency')
        end
    end
    
    
    % Plot TFF by stimulus frequency
    figure(5)
    subplot(2,1,1)
    hold on
    if x==3
         errorbar(A([1:3 5]),ttf(x).ttf_FrM([1:3 5],:),ttf(x).ttf_FrM([1:3 5],:)-ttf(x).ttf_FrCI([1:3 5],1),ttf(x).ttf_FrCI([1:3 5],2)-ttf(x).ttf_FrM([1:3 5],:),['-o' color])
    else
         errorbar(A,ttf(x).ttf_FrM,ttf(x).ttf_FrM-ttf(x).ttf_FrCI(:,1),ttf(x).ttf_FrCI(:,2)-ttf(x).ttf_FrM,['-o' color])
    end
   
    title(observerID)
    ylabel('power spectra for stimulus frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[0 0.02];
  
    % plot background across channels
    if x==3
        f=Fs*(0:((Fs*dur_in_sec)/2))/(Fs*dur_in_sec);
        background=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        background2=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        counter=0;
        for yyy=1:size(background,1)-1
            if isnan(background(yyy))==1
                background2=cat(1,background2(1:yyy-1-counter,:),background(yyy+1:end,:));
                counter=counter+1;
            end
        end
        rand_trial=sort(randi(size(background2,1),1,21));
        background2=background2(rand_trial,:);
        backgroundM=nanmean(background2,1);
        ft=fft(backgroundM);
        P=abs(ft/(Fs*dur_in_sec));
        ttf_BKGD=P(1:(Fs*dur_in_sec)/2+1);
        Bootstat=bootstrp(1000,@nanmean,background,1);
        for yy=1:size(Bootstat,1)
            boot_ft=fft(Bootstat(yy,:));
            P_boot=abs(boot_ft/(Fs*dur_in_sec));
            ttf_bkgd_boot(yy,:)=P_boot(1:(Fs*dur_in_sec)/2+1);
        end
        
        ttf_bkgd_boot=sort(ttf_bkgd_boot,1);
        ttf_bkgdCI=ttf_bkgd_boot([50 950],:);
        
        for bb=1:length(A)
            temp=find(ttf(x).f>=A(bb));
            ttf_bkgd_Fr(bb,:)=ttf_BKGD(:,temp(1));
            ttf_bkgdCI_Fr(bb,:)=ttf_bkgdCI(:,temp(1));
        end
        
        temp2=find(ttf(x).f>=0.5 & ttf(x).f<=4);
        ttf_bkgd_Fr(1,:)=ttf_BKGD(:,temp(1));
        ttf_bkgdCI_Fr(1,:)=ttf_bkgdCI(:,temp(1));
     
        errorbar(A,ttf_bkgd_Fr,ttf_bkgdCI_Fr(:,1),ttf_bkgdCI_Fr(:,2),'-o','Color',[0.5 0.5 0.5])
    end
    
    % Plot TFF for delta frequency
    figure(22)
    subplot(1,2,1)
    hold on
    if x==3
         errorbar(A([1:3 5]),ttf(x).ttf_deltaM([1:3 5],:),ttf(x).ttf_deltaM([1:3 5],:)-ttf(x).ttf_deltaCI([1:3 5],1),ttf(x).ttf_deltaCI([1:3 5],2)-ttf(x).ttf_deltaM([1:3 5],:),['-o' color])
    else
         errorbar(A,ttf(x).ttf_deltaM,ttf(x).ttf_deltaM-ttf(x).ttf_deltaCI(:,1),ttf(x).ttf_deltaCI(:,2)-ttf(x).ttf_deltaM,['-o' color])
    end
    
   
    title('Delta frequency (0.5-4Hz)')
    xlabel('temporal frequency of stimulus')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[0 0.02];
    
    
    % Plot TFF for alpha frequency
    figure(22)
    subplot(1,2,2)
    hold on
    if x==3
         errorbar(A([1:3 5]),ttf(x).ttf_alphaM([1:3 5],:),ttf(x).ttf_alphaM([1:3 5],:)-ttf(x).ttf_alphaCI([1:3 5],1),ttf(x).ttf_alphaCI([1:3 5],2)-ttf(x).ttf_alphaM([1:3 5],:),['-o' color])
    else
         errorbar(A,ttf(x).ttf_alphaM,ttf(x).ttf_alphaM-ttf(x).ttf_alphaCI(:,1),ttf(x).ttf_alphaCI(:,2)-ttf(x).ttf_alphaM,['-o' color])
    end
   
    title(observerID)
    title('Alpha frequency (8-12Hz)')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[0 0.02];
    
    
    % Plot superimposed luminance, red/green, and blue/yellow in time
    % domain
    figure(11)
    for z=1:length(A)
        subplot(3,2,z)
        plot(XX,squeeze(nanmean(processedVEPdata(x).vep_Fr(z,:,:),2)),['-' color])
        title(['frequency=' num2str(A(z))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 parsedVEPdata(x).dur_in_freq/parsedVEPdata(x).Fs];
        hold on
    end
    
    
    
%     % info to save per channel
%     VEP_main.mtrp.group;
%     VEP_main.mtrp.observer;

end



% save(filenameComp)
%end

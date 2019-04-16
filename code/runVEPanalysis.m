
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
nulling(1)=input('LM nulling value:');
nulling(2)=input('S nulling value:');

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
    
    [parsedVEPdata(x)]=parseVEP(VEP_main,'dur_in_sec',dur_in_sec,'starttime',starttime,'bandstop60',true,'plot_sessions',false);
    Fs=VEP_main(1).vepDataStruct.params.frequencyInHz;
    XX=(1:length(parsedVEPdata(x).vep_Fr))/Fs;
    A=unique(VEP_main(x).mtrp.TFtrials);
end

%% Process VEP data (gets rid of poor quality trials, and normalizes signal)
    [processedVEPdata]=preprocessVEP(parsedVEPdata,'dur_in_sec',dur_in_sec,'plot_all',false);
    
    % analyze background data across channels
    [ttf_bkgd]=vepBKGD(processedVEPdata,Fs,dur_in_sec,A);
    
    [fooof_bkgd,fooof_bkgd5,fooof_bkgd95]=runFOOOF(ttf_bkgd.ttf_BKGD,ttf_bkgd.ttf_bkgdCI,Fs,dur_in_sec);
%%
for x=1:3
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
    %% Calculate TTF
    [ttf(x)]=calcVEPttf(processedVEPdata(x).vep_Fr,'dur_in_sec',dur_in_sec,'plot_all',false,'TemporalFrequency',A);

    %% FOOOF
    [fooof_results(x,:),fooof_results5(x,:),fooof_results95(x,:)]=runFOOOF(ttf(x).ttf_M,ttf(x).ttf_CI,Fs,dur_in_sec);

    
    %% Plotting
    % Plot median Visual discomfort data
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
    ax.YLim=[0 11];
    
    % Plot TTF
    if x==1
        figTTF=figure('Name',observerID);
    end
    figure(figTTF)
    hold on
    for YY=1:length(A)
        YYY=x+((YY-1)*3);
        subplot(5,3,YYY)
        xdata=cat(2,ttf(x).f,fliplr(ttf(x).f));
        ydata=cat(2,squeeze(ttf(x).ttf_CI(YY,:,1)),fliplr(squeeze(ttf(x).ttf_CI(YY,:,2))));
        TEMP=fill(xdata,ydata,Color,'EdgeColor','none');
        hold on
        plot(ttf(x).f,ttf(x).ttf_M(YY,:),['-' color])
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YScale='log';
        ax.XLim=[min(ttf(x).f) 100];
        ax.XTick=A;
        ax.YLim=[0.00001 0.2];
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
        neg=ttf(x).ttf_FrM([1:3 5],:)-ttf(x).ttf_FrCI([1:3 5],1);
        pos=ttf(x).ttf_FrCI([1:3 5],2)-ttf(x).ttf_FrM([1:3 5],:);
        errorbar(A([1:3 5]),ttf(x).ttf_FrM([1:3 5],:),neg,pos,['-o' color])
        
        neg=ttf_bkgd.ttf_bkgd_Fr-ttf_bkgd.ttf_bkgdCI_Fr(:,1);
        pos=ttf_bkgd.ttf_bkgdCI_Fr(:,2)-ttf_bkgd.ttf_bkgd_Fr;
        errorbar(A,ttf_bkgd.ttf_bkgd_Fr,neg,pos,'-o','Color',[0.5 0.5 0.5])
    else
        neg=ttf(x).ttf_FrM(:,:)-ttf(x).ttf_FrCI(:,1);
        pos=ttf(x).ttf_FrCI(:,2)-ttf(x).ttf_FrM(:,:);
        errorbar(A,ttf(x).ttf_FrM,neg,pos,['-o' color])
    end
   
    title(observerID)
    ylabel('power spectra for stimulus frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[-0.002 0.02];
  
    % plot background across channels
    if x==3     
        errorbar(A,ttf_bkgd.ttf_bkgd_Fr,ttf_bkgd.ttf_bkgdCI_Fr(:,1),ttf_bkgd.ttf_bkgdCI_Fr(:,2),'-o','Color',[0.5 0.5 0.5])
    end
   
     % Plot psd FOOOF
    figure(7)
     subplot(1,2,1)
    hold on
    if x==3
        neg=ttf(x).ttf_FrM([1:3 5],:)-ttf(x).ttf_FrCI([1:3 5],1);
        pos=ttf(x).ttf_FrCI([1:3 5],2)-ttf(x).ttf_FrM([1:3 5],:);
        errorbar(A([1:3 5]),ttf(x).ttf_FrM([1:3 5],:),neg,pos,['-o' color])
        
        neg=ttf_bkgd.ttf_bkgd_Fr-ttf_bkgd.ttf_bkgdCI_Fr(:,1);
        pos=ttf_bkgd.ttf_bkgdCI_Fr(:,2)-ttf_bkgd.ttf_bkgd_Fr;
        errorbar(A,ttf_bkgd.ttf_bkgd_Fr,neg,pos,'-o','Color',[0.5 0.5 0.5])
    else
        neg=ttf(x).ttf_FrM(:,:)-ttf(x).ttf_FrCI(:,1);
        pos=ttf(x).ttf_FrCI(:,2)-ttf(x).ttf_FrM(:,:);
        errorbar(A,ttf(x).ttf_FrM,neg,pos,['-o' color])
    end
    
    title(observerID)
    ylabel('power spectra for stimulus frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[-0.002 0.02];
    
    
    % get FOOOF peak psd 
    
    for a=1:length(A)
        xdata=fooof_results(x,a).freqs;
        ydata=10.^(fooof_results(x,a).power_spectrum)-10.^(fooof_results(x,a).bg_fit);
        ydata5=10.^(fooof_results5(x,a).power_spectrum)-10.^(fooof_results5(x,a).bg_fit);
        ydata95=10.^(fooof_results95(x,a).power_spectrum)-10.^(fooof_results95(x,a).bg_fit);
       
        peak_freq=A(a);
        temp=abs(xdata-peak_freq);
        temp2=find(temp==min(temp));
        if length(temp2)>1
            temp3=max(ydata(temp2));
            peak_freq_loc(a)=find(ydata==temp3);
        else
            peak_freq_loc(a)=temp2;
        end
        
        
        fooof_peak_Fr(x,a)=ydata(peak_freq_loc(a));
        fooof_peak_Fr5(x,a)=ydata5(peak_freq_loc(a));
        fooof_peak_Fr95(x,a)=ydata95(peak_freq_loc(a));
        
    end
    
    % get harmonics
    harmonic=[1 2 3 4];
    for a=1:length(A)
        xdata=fooof_results(x,a).freqs;
        ydata=10.^(fooof_results(x,a).power_spectrum)-10.^(fooof_results(x,a).bg_fit);
        
        peak_freq=A(a);
        peak_freq_harmonic=harmonic.*peak_freq;
        peak_freq_harmonic=peak_freq_harmonic(find(peak_freq_harmonic<100 & peak_freq_harmonic>1));
        for b=1:length(peak_freq_harmonic)
            temp=abs(xdata-peak_freq_harmonic(b));
            temp2=find(temp==min(temp));
            if length(temp2)>1
                temp3=max(ydata(temp2));
                peak_freq_harm_loc(b)=find(ydata==temp3);
            else
                peak_freq_harm_loc(b)=temp2;
            end
        end
        figure(12)
        plot(xdata,ydata,'k')
        hold on
        plot(xdata(peak_freq_harm_loc),ydata(peak_freq_harm_loc),'or')
        ax=gca;
        ax.Box='off';
        ax.TickDir='out';
        ax.XLim=[0 100];
        ax.YLim=[-0.002 0.02];
        pause
        hold off
        
        fooof_peak_harmonics{x,a,:}=ydata(peak_freq_harm_loc);
        fooof_peak_harmonics_freq{x,a,:}=xdata(peak_freq_harm_loc);
        
        clear ydata xdata peak_freq_harm_loc peak_freq_harmonic peak_freq
    end
    
    
        
    
    if x==3
        xdata=fooof_bkgd.freqs;
        ydata=10.^(fooof_bkgd.power_spectrum)-10.^(fooof_bkgd.bg_fit);
        ydata5=10.^(fooof_bkgd5.power_spectrum)-10.^(fooof_bkgd5.bg_fit);
        ydata95=10.^(fooof_bkgd95.power_spectrum)-10.^(fooof_bkgd95.bg_fit);

        for a=1:length(A)
           fooof_bkgdFr(:,a)=ydata(peak_freq_loc(a));
           fooof_bkgdFr5(:,a)=ydata5(peak_freq_loc(a));
           fooof_bkgdFr95(:,a)=ydata95(peak_freq_loc(a));
        end
    end
    
    figure(7)
    subplot(1,2,2)
    hold on
    if x==3
        errorbar(A([1:3 5]),fooof_peak_Fr(x,[1:3 5]),fooof_peak_Fr(x,[1:3 5])-fooof_peak_Fr5(x,[1:3 5]),fooof_peak_Fr95(x,[1:3 5])-fooof_peak_Fr(x,[1:3 5]),['-o' color])

        errorbar(A,fooof_bkgdFr,fooof_bkgdFr-fooof_bkgdFr5,fooof_bkgdFr95-fooof_bkgdFr,'-o','Color',[0.5 0.5 0.5])
         
    else
        errorbar(A,fooof_peak_Fr(x,:),fooof_peak_Fr(x,:)-fooof_peak_Fr5(x,:),fooof_peak_Fr95(x,:)-fooof_peak_Fr(x,:),['-o' color])
    end
   
    title([observerID ' fooofed'])
    ylabel('power spectra for stimulus frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[-0.002 0.02];
  
    
    % Plot superimposed luminance, red/green, and blue/yellow in time
    % domain
    figure(11)
    for z=1:length(A)
        subplot(3,2,z)
        vep_temp=squeeze(nanmedian(processedVEPdata(x).vep_Fr(z,:,:),2));
        plot(XX,vep_temp,['-' color])
        title(['frequency=' num2str(A(z))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 parsedVEPdata(x).dur_in_freq/parsedVEPdata(x).Fs];
        hold on
    end

    vds(x,:,:)=parsedVEPdata(x).vds_Fr;
    vep_Fr(x,:,:,:)=processedVEPdata(x).vep_Fr;
    vep_BKGD(x,:,:)=ttf_bkgd.vep_bkgd;
end

compiledData.observerID=observerID;
compiledData.group=VEP_main(1).mtrp.group;
compiledData.Fs=Fs;
compiledData.vds=vds;
compiledData.vep_Fr=vep_Fr;
compiledData.vep_bkgd=vep_BKGD;
compiledData.fooof_peak_Fr=fooof_peak_Fr;
compiledData.fooof_bkgd_Fr=fooof_bkgdFr;
compiledData.fooof_peak_harmonics=fooof_peak_harmonics;
compiledData.fooof_peak_harmonics_freq=fooof_peak_harmonics_freq;
compiledData.fooof_results=fooof_results;
compiledData.fooof_bkgd=fooof_bkgd;
compiledData.ttf_M=ttf(x).ttf_M;
compiledData.ttf_CI=ttf(x).ttf_CI;
compiledData.ttf_bkgd=ttf_bkgd;
compiledData.nulling=nulling;

save(filenameComp,'compiledData')


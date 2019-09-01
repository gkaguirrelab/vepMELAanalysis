
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
                Color=[0 0 0];
            case 2
                expID='LM004';
                color='r';
                Color=[1 0 0];
            case 3
                expID='S0004';
                color='b';
                Color=[0 0 1];
    end
    
    filenameMAT=fullfile(getpref('vepMELAanalysis','melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID],['Exp' expID '_' observerID 'compiled.mat']);

    open(filenameMAT);
    VEP_main=ans.VEP;

    %% Parse VEP data
    
    [parsedVEPdata(x)]=parseVEP(VEP_main,'dur_in_sec',dur_in_sec,'starttime',starttime,'bandstop60',true,'plot_sessions',false);
    Fs=VEP_main(1).vepDataStruct.params.frequencyInHz;
    XX=(1:length(parsedVEPdata(x).vep_Fr))/Fs;
    TemporalFrequency=unique(VEP_main(x).mtrp.TFtrials);
end

%% Process VEP data (gets rid of poor quality trials, and normalizes signal)
    [processedVEPdata]=preprocessVEP(parsedVEPdata,'dur_in_sec',dur_in_sec,'plot_all',false);

%%
for x=1:3
    switch x
            case 1
                expID='LMS05';
                color='k';
                Color=[0 0 0];
                x0=[1 2 1];
                flicker_freq=1:5;
            case 2
                expID='LM004';
                color='r';
                Color=[1 0 0];
                 x0=[2 2 1];
                flicker_freq=1:5;
            case 3
                expID='S0004';
                color='b';
                Color=[0 0 1];
                x0=[6 2 1];
                flicker_freq=[1:3 5];
    end
    %% Calculate TTF
    [ttf(x)]=calcVEPttf(processedVEPdata(x).vep_Fr,'dur_in_sec',dur_in_sec,'plot_all',false,'TemporalFrequency',TemporalFrequency);

    %% FOOOF
    [fooof_results(x,:),fooof_results5(x,:),fooof_results95(x,:)]=runFOOOF(ttf(x).ttf_M,ttf(x).ttf_CI,'plot_all',true);

    
    %% Plotting
    % Plot median Visual discomfort data
    figure(5)
    VDSm=nanmedian(parsedVEPdata(x).vds_Fr,2);
    VDSstd=nanstd(parsedVEPdata(x).vds_Fr,[],2);
    errorbar(TemporalFrequency,VDSm,VDSstd,['-o' color])
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
    for YY=1:length(TemporalFrequency)
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
        ax.XTick=TemporalFrequency;
        ax.YLim=[0.00001 0.2];
        title(num2str(TemporalFrequency(YY)))
        if YYY==7
            ylabel('Power')
        end
        
        if YYY==14
            xlabel('Frequency')
        end
    end
    
   
     % Plot psd FOOOF
    figure(7)
    subplot(1,2,1)
    hold on
    markerline=['-o' color];markeredge=Color;markerface=Color;
    plotWithErrorbars(TemporalFrequency(flicker_freq)',ttf(x).ttf_FrM(flicker_freq)',ttf(x).ttf_FrCI(flicker_freq,:)',markerline,markeredge,markerface)
    title(observerID)
    ylabel('amplitude of stimulus frequency (mV)')
    ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.04];
    
    
    % get FOOOF peak psd 
    sbplot1=1:2:10; sbplot2=2:2:10;
    
    for a=1:length(TemporalFrequency)
        xdata=fooof_results(x,a).freqs;
        ydata_psd=fooof_results(x,a).power_spectrum;
        ydata_ap=fooof_results(x,a).bg_fit;
        ydata=10.^(fooof_results(x,a).power_spectrum)-10.^(fooof_results(x,a).bg_fit);
        ydata5=10.^(fooof_results5(x,a).power_spectrum)-10.^(fooof_results5(x,a).bg_fit);
        ydata95=10.^(fooof_results95(x,a).power_spectrum)-10.^(fooof_results95(x,a).bg_fit);

        peak_freq=TemporalFrequency(a);
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
        
        figure(15)
        subplot(5,2,sbplot1(a))
        plot(xdata,10.^ydata_psd,'k')
        hold on
        plot(xdata,10.^ydata_ap,'--g')
        ylabel('amplitude (mV)')
        xlabel('frequency (Hz)')
        ax=gca;ax.XLim=[0 60];ax.YLim=[-0.002 0.04];ax.Box='off';ax.TickDir='out';
       
        subplot(5,2,sbplot2(a))
        plot(xdata,ydata,'k')
        ylabel('amplitude (mV)')
        xlabel('frequency (Hz)')
        ax=gca;ax.XLim=[0 60];ax.YLim=[-0.002 0.04];ax.Box='off';ax.TickDir='out';
        
    end
    pause
    
    % get harmonics
    harmonic={[1.625 3.25 4.875 6.5 8.125];[3.25 6.5 9.75 13 16.25];[7.5 15 22.5 30 37.5];[15 30 45 75 90];[30 90]};
    for a=1:length(TemporalFrequency)
        xdata=fooof_results(x,a).freqs;
        ydata=10.^(fooof_results(x,a).power_spectrum)-10.^(fooof_results(x,a).bg_fit);
        
        peak_freq_harmonic=cell2mat(harmonic(a,:));
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
        ax.YLim=[-0.002 0.04];
%         pause
        hold off
        
        fooof_peak_harmonics{x,a,:}=ydata(peak_freq_harm_loc);
        fooof_peak_harmonics_freq{x,a,:}=xdata(peak_freq_harm_loc);
        
        clear ydata xdata peak_freq_harm_loc peak_freq_harmonic peak_freq
    end
         
    figure(7)
    subplot(1,2,2)
    hold on
    markerline=['-o' color];markeredge=Color;markerface=Color;
    plotWithErrorbars(TemporalFrequency(flicker_freq),fooof_peak_Fr(x,flicker_freq),cat(1,fooof_peak_Fr5(x,flicker_freq),fooof_peak_Fr95(x,flicker_freq)),markerline,markeredge,markerface)
    title('FOOOF')
    ylabel('amplitude of stimulus frequency (mV)')
    ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.04];


    
    % Plot superimposed luminance, red/green, and blue/yellow in time
    % domain
    figure(11)
    for z=1:length(TemporalFrequency)
        subplot(3,2,z)
        vep_temp=squeeze(nanmedian(processedVEPdata(x).vep_Fr(z,:,:),2));
        plot(XX,vep_temp,['-' color])
        title(['frequency=' num2str(TemporalFrequency(z))]);
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
end

compiledData.observerID=observerID;
compiledData.group=VEP_main(1).mtrp.group;
compiledData.Fs=Fs;
compiledData.vds=vds;
compiledData.vep_Fr=vep_Fr;
compiledData.fooof_peak_Fr=fooof_peak_Fr;
compiledData.fooof_peak_harmonics=fooof_peak_harmonics;
compiledData.fooof_peak_harmonics_freq=fooof_peak_harmonics_freq;
compiledData.fooof_results=fooof_results;
compiledData.ttf_M=ttf(x).ttf_M;
compiledData.ttf_CI=ttf(x).ttf_CI;
compiledData.nulling=nulling;

save(filenameComp,'compiledData')


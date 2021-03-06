% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169';...
    'MELA_0194';'MELA_0179';...
    'MELA_0191';'MELA_0174';...
    'MELA_0120';'MELA_0171';...
    'MELA_0201';'MELA_0207';...
    'MELA_0209'];

counter_MVA=0;
counter_HAF=0;
A=[1.625 3.25 7.5 15 30];
lb=50;
ub=950;

savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

for x=1:size(subjects,1)
    
        observerID=subjects(x,:);
        savePath=fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
        'vepMELAanalysis','allChannels');

        filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);

        open(filenameComp);
        compiledData=ans.compiledData;

        if compiledData(1).group=='MWVA'
            counter_MVA=counter_MVA+1;
            MVA_ID(counter_MVA,:)=compiledData.observerID;
            MVA_vds(counter_MVA,:,:,:)=compiledData.vds;
            MVA_fooof_fr(counter_MVA,:,:)=compiledData.fooof_peak_Fr;
            MVA_nulling(counter_MVA,:)=compiledData.nulling;
            for xx=1:3
                for yy=1:5
                    MVA_ap(counter_MVA,xx,yy,:)=compiledData.fooof_results(xx,yy).bg_fit;
                    MVA_psd(counter_MVA,xx,yy,:)=compiledData.fooof_results(xx,yy).power_spectrum;
                    MVA_fooof(counter_MVA,xx,yy,:)=compiledData.fooof_results(xx,yy).fooofed_spectrum;
                end
            end
        else if compiledData(1).group=='HA f'
                counter_HAF=counter_HAF+1;
                HAF_ID(counter_HAF,:)=compiledData.observerID;
                HAF_vds(counter_HAF,:,:,:)=compiledData.vds;
                HAF_fooof_fr(counter_HAF,:,:)=compiledData.fooof_peak_Fr;
                HAF_nulling(counter_HAF,:)=compiledData.nulling;
                   for xx=1:3
                        for yy=1:5
                            HAF_ap(counter_HAF,xx,yy,:)=compiledData.fooof_results(xx,yy).bg_fit;
                            HAF_psd(counter_HAF,xx,yy,:)=compiledData.fooof_results(xx,yy).power_spectrum;
                            HAF_fooof(counter_HAF,xx,yy,:)=compiledData.fooof_results(xx,yy).fooofed_spectrum;
                        end
                    end
            else
                disp('error: group');
            end
        end

end

freqs=compiledData.fooof_results(1,1).freqs;
clear x xx yy ans


% Calculate bootstrapped confidence intervals for fooof by frequency
for a=1:size(HAF_vds,2)
    for b=1:size(HAF_vds,3)
        MVA_temp=squeeze(MVA_fooof_fr(:,a,b));
        Bootstat=bootstrp(1000,@nanmedian,MVA_temp,1);
        Bootstat=sort(Bootstat,1);
        MVA_fooof_CI(a,b,:)=Bootstat([lb ub],:)';
        
        HAF_temp=squeeze(HAF_fooof_fr(:,a,b));
        Bootstat2=bootstrp(1000,@nanmedian,HAF_temp,1);
        Bootstat2=sort(Bootstat2,1);
        HAF_fooof_CI(a,b,:)=Bootstat2([lb ub],:)';
        
        MVA_temp2=squeeze(squeeze(squeeze(nanmedian(MVA_vds(:,a,b,:),4))));
        Bootstat3=bootstrp(1000,@nanmedian,MVA_temp2,1);
        Bootstat3=sort(Bootstat3,1);
        MVA_vds_CI(a,b,:)=Bootstat3([lb ub],:)';
        
        HAF_temp2=squeeze(squeeze(squeeze(nanmedian(HAF_vds(:,a,b,:),4))));
        Bootstat4=bootstrp(1000,@nanmedian,HAF_temp2,1);
        Bootstat4=sort(Bootstat4,1);
        HAF_vds_CI(a,b,:)=Bootstat4([lb ub],:)';
    end
end

clear a b
        
    % Plot TFF by stimulus frequency
    figure(5)
    MVA_fooof_M=squeeze(nanmedian(MVA_fooof_fr,1));
    HAF_fooof_M=squeeze(nanmedian(HAF_fooof_fr,1));
    for y=1:size(HAF_vds,2)
        switch y
            case 1
                color='k';
                flicker_stim=1:5;
            case 2
                color='r';
                flicker_stim=1:5;
            case 3
                color='b';
                flicker_stim=[1:3 5];
        end
        
        subplot(2,2,1)
        neg=MVA_fooof_M-squeeze(MVA_fooof_CI(:,:,1));
        pos=squeeze(MVA_fooof_CI(:,:,2))-MVA_fooof_M;
        errorbar(A(flicker_stim),MVA_fooof_M(y,flicker_stim),neg(y,flicker_stim),pos(y,flicker_stim),['-o' color],'LineWidth',2,'MarkerFaceColor',color)
        hold on
        title(['Migraine with visual aura (n=' num2str(counter_MVA) ')'])
        ylabel('amplitude at stimulus frequency (mV)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[-0.001 0.03];
        
        subplot(2,2,2)
        neg=HAF_fooof_M-squeeze(HAF_fooof_CI(:,:,1));
        pos=squeeze(HAF_fooof_CI(:,:,2))-HAF_fooof_M;
        errorbar(A(flicker_stim),HAF_fooof_M(y,flicker_stim),neg(y,flicker_stim),pos(y,flicker_stim),['-o' color],'LineWidth',2)
        hold on
        title(['Headache free (n=' num2str(counter_HAF) ')'])
        ylabel('amplitude at stimulus frequency (mV)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[-0.001 0.03];
        

        
    end



   % plot individual subjects fooof frequency
for y=1:size(HAF_vds,2)
        switch y
            case 1
                color='k';
                flicker_stim=1:5;
            case 2
                color='r';
                flicker_stim=1:5;
            case 3
                color='b';
                flicker_stim=[1:3 5];
        end
            figure(101)
            for c=1:size(HAF_fooof_fr,1)
                subplot(3,3,c)
                hold on
                plot(A(flicker_stim),squeeze(HAF_fooof_fr(c,y,flicker_stim)),['-o' color],'MarkerFaceColor',color)
                title(num2str(HAF_ID(c,:)))
                ax=gca;
                ax.TickDir='out';
                ax.Box='off';
                ax.XScale='log';
                ax.XLim=[0.95 35];
                ax.YLim=[-0.001 0.03];
            end
            
            figure(100)
            for b=1:size(MVA_fooof_fr,1)
                subplot(3,3,b)
                hold on
                plot(A(flicker_stim),squeeze(MVA_fooof_fr(b,y,flicker_stim)),['-o' color],'MarkerFaceColor',color)
                title(num2str(MVA_ID(b,:)))
                ylabel('amplitude at stimulus frequency (mV)')
                ax=gca;
                ax.TickDir='out';
                ax.Box='off';
                ax.XScale='log';
                ax.XLim=[0.95 35];
                ax.YLim=[-0.001 0.03];
            end
end
    
    % Plot median Visual discomfort data
    figure(5)
    MVA_vds=squeeze(nanmedian(MVA_vds,4));
    HAF_vds=squeeze(nanmedian(HAF_vds,4));
    for y=1:size(HAF_vds,2)
        switch y
            case 1
                color='k';
                flicker_stim=1:5;
            case 2
                color='r';
                flicker_stim=1:5;
            case 3
                color='b';
                flicker_stim=[1:3 5];
        end
        subplot(2,2,3)
        VDS_MVA=squeeze(nanmedian(MVA_vds,1));
        neg=VDS_MVA-squeeze(MVA_vds_CI(:,:,1));
        pos=squeeze(MVA_vds_CI(:,:,2))-VDS_MVA;
        errorbar(A(flicker_stim),VDS_MVA(y,flicker_stim),neg(y,flicker_stim),pos(y,flicker_stim),['-o' color],'MarkerFaceColor',color,'LineWidth',2)
        hold on
        ylabel('flicker discomfort')
        xlabel('temporal frequency of stimulus')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[0 11];
    
        subplot(2,2,4)
        VDS_HAF=squeeze(nanmedian(HAF_vds,1));
        neg=VDS_HAF-squeeze(HAF_vds_CI(:,:,1));
        pos=squeeze(HAF_vds_CI(:,:,2))-VDS_HAF;
        errorbar(A(flicker_stim),VDS_HAF(y,flicker_stim),neg(y,flicker_stim),pos(y,flicker_stim),['-o' color],'MarkerFaceColor','w','LineWidth',2)
        hold on
        ylabel('flicker discomfort')
        xlabel('temporal frequency of stimulus')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[0 11];
    end

    % Plot aperiodic signal
    figure(9)
    for xx=1:size(MVA_ap,1)
        MVA_ap_subject(xx,:)=nanmedian(reshape(squeeze(MVA_ap(xx,:,:,:)),15,size(MVA_ap,4)));
        
         % plot median by subject
            subplot(1,3,1)
            hold on
            plot(freqs,MVA_ap_subject(xx,:),'k')
            title('Migraine with visual aura')
            ylabel('amplitude at stimulus frequency (mV)')
            ax=gca;
            ax.TickDir='out';
            ax.Box='off';
            ax.XLim=[0 35];
            ax.YLim=[-4 -2];
    end
    
    for xx=1:size(HAF_ap,1)
        HAF_ap_subject(xx,:)=nanmedian(reshape(squeeze(HAF_ap(xx,:,:,:)),15,size(HAF_ap,4)));
        
         % plot median by subject
            subplot(1,3,2)
            hold on
            plot(freqs,HAF_ap_subject(xx,:),'--k')
            title('Headache free')
            ylabel('amplitude at stimulus frequency (mV)')
            ax=gca;
            ax.TickDir='out';
            ax.Box='off';
            ax.XLim=[0 35];
            ax.YLim=[-4 -2];
    end
    
   
    
    subplot(1,3,3)
    hold on
    
    Bootstat=bootstrp(1000,@nanmedian,HAF_ap_subject,1);
    Bootstat=sort(Bootstat,1);
    Bootstat=Bootstat([lb ub],:);
    median=nanmedian(HAF_ap_subject,1);
    Y=cat(2,Bootstat(2,:),fliplr(Bootstat(1,:)));
    X=cat(2,freqs,fliplr(freqs));
    TEMP=fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    plot(freqs,median,'--k')
    
    Bootstat=bootstrp(1000,@nanmedian,MVA_ap_subject,1);
    Bootstat=sort(Bootstat,1);
    Bootstat=Bootstat([lb ub],:);
    median=nanmedian(MVA_ap_subject,1);
    Y=cat(2,Bootstat(2,:),fliplr(Bootstat(1,:)));
    X=cat(2,freqs,fliplr(freqs));
    TEMP=fill(X,Y,[0.5 0.5 0.5],'EdgeColor','none');
    plot(freqs,median,'k')
    
    title('Comparison')
    ylabel('amplitude at stimulus frequency (mV)')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 35];
    ax.YLim=[-4 -2];
            
  
    
    % Plot nulling
    figure(10)
    plot([1 2],MVA_nulling,'ok','MarkerFaceColor','k')
    hold on
    plot([1 2],HAF_nulling,'ok')
    plot([0 1 3],[0 0 0],'--')
    title('Nulling values')
    ylabel('RGB adjustment values')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XTick=[1 2];
    ax.XTickLabel={'L-M','S'};
    ax.XLim=[0.5 2.5];
    ax.YLim=[-0.12 0.02];
    
% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169';...
    'MELA_0194';'MELA_0179';...
    'MELA_0191';'MELA_0174';...
    'MELA_0120';'MELA_0171';...
    'MELA_0201';'MELA_0204';...
    'MELA_0207';'MELA_0209';...
    'MELA_0213';'MELA_0208'];


counter_MVA=0;
counter_HAF=0;
TemporalFrequency=[1.625 3.25 7.5 15 30];
lb=50;
ub=950;

savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

load('/Users/carlynpattersongentile/Documents/MATLAB/Carlyn/vepMELA_subjectInfo.mat')

for x=1:size(subjects,1)
    
        observerID=subjects(x,:);
        savePath=fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
        'vepMELAanalysis','allChannels');

        filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);

        open(filenameComp);
        compiledData=ans.compiledData;
        temp=find(strcmp(table2array(scoreTable(:,1)),observerID));
        temp2=scoreTable(temp,:);
        compiledData.subject=temp2;
        compiledData_all(x,:)=compiledData;
        
        if compiledData(1).group=='MWVA'
            counter_MVA=counter_MVA+1;
            compiledData_MVA(counter_MVA,:)=compiledData;
            scoreTable_MVA(counter_MVA,:)=temp2;
                else if compiledData(1).group=='HA f'
                    counter_HAF=counter_HAF+1;
                   compiledData_HAF(counter_HAF,:)=compiledData;
                   scoreTable_HAF(counter_HAF,:)=temp2;
                    else
                    disp('error: group');
                end
        end
end

clear ans compiledData counter_HAF counter_MVA observerID x filenameComp

% flicker discomfort
VDS=calcVDS(compiledData_all,lb,ub);
VDS_mva=calcVDS(compiledData_MVA,lb,ub);
VDS_haf=calcVDS(compiledData_HAF,lb,ub);

% Plot the harmonics across subjects
[pred60lms, pred60lm, pred60s]=plotFooofHarmonics(compiledData_all,3,lb,ub);

% Plot the sum of harmonics across subjects
[LMS, LM, S, LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_all,lb,ub, pred60lms, pred60lm, pred60s);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[2 1 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[2 3 0.5]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(1)
subplot(3,2,1)
hold on
errorbar(TemporalFrequency,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Sum of harmonics'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.035];

subplot(3,2,3)
hold on
errorbar(TemporalFrequency,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.035];

subplot(3,2,5)
hold on
errorbar(TemporalFrequency([1:3 5]),Sm([1:3 5]),Sm([1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm([1:3 5]),'ob','LineWidth',2)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.035];

% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
subplot(3,2,[2 4 6])
hold on
plot(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMSM,VDS.LMM,VDS.SM([1:3 5])),'.','Color',[0.5 0.5 0.5])
lsline
plot(LMSm,VDS.LMSM,'ok','MarkerFaceColor','k','MarkerSize',8)
plot(LMm,VDS.LMM,'or','MarkerFaceColor','r','MarkerSize',8)
plot(Sm([1:3 5]),VDS.SM([1:3 5]),'ob','MarkerFaceColor','b','MarkerSize',8)
[Rmva_vdsvep,Pmva_vdsvep]=corrcoef(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMSM,VDS.LMM,VDS.SM([1:3 5])));
title(['n=14, R squared=' num2str(Rmva_vdsvep(1,2)^2) ' , p=' num2str(Pmva_vdsvep(1,2))])
ylabel('flicker discomfort')
xlabel('amplitude at stimulus frequency (mV)')
ax=gca; ax.XLim=[-0.002 0.03]; ax.YLim=[0 7]; ax.TickDir='out'; ax.Box='off';

% calculate ANOVA for flicker discomfort by group, post-receptoral pathway,
% and temporal frequency
VEP_all=[];
VDS_all=[];
group=[];
PRP=[];
TF=[];

[pred60lms, pred60lm, pred60s]=plotFooofHarmonics(compiledData_MVA,3,lb,ub);
[LMS, LM, S, LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_MVA,lb,ub, pred60lms, pred60lm, pred60s);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 1 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(2)
subplot(3,2,1)
hold on
errorbar(TemporalFrequency,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Sum of harmonics'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

subplot(3,2,3)
hold on
errorbar(TemporalFrequency,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

subplot(3,2,5)
hold on
errorbar(TemporalFrequency([1:3 5]),Sm([1:3 5]),Sm([1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm([1:3 5]),'ob','LineWidth',2)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

VEP_harmonic=cat(3,LMS,LM,S);
for i=1:size(compiledData_MVA,1)
    temp=nanmedian(compiledData_MVA(i).vds,3);
    temp2=squeeze(VEP_harmonic(i,:,:))';
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            VDS_all=cat(2,VDS_all,temp(j,k));
            VEP_all=cat(2,VEP_all,temp2(j,k));
            group=cat(2,group,{'MVA'});
            TF=cat(2,TF,TemporalFrequency(k));
        
            switch j
                case 1
                    PRP=cat(2,PRP,{'LMS'});
                case 2
                    PRP=cat(2,PRP,{'LM'});
                case 3
                    PRP=cat(2,PRP,{'S'});
            end
        end
    end
end

clear temp temp2

% plot luminance 30 Hz VEP response as a function of headache frequency
figure(6)
for i=1:size(compiledData_MVA,1)
    X_mva(i)=table2array(compiledData_MVA(i).subject(:,6));
    Y_mva(i)=LMS(i,5);
end

[Rmva_vepHAf,Pmva_vepHAf]=corrcoef(X_mva,Y_mva);
plot(X_mva,Y_mva,'ok','MarkerFaceColor','k','MarkerSize',8)
hold on
lsline
markerline='ok';markeredge=[0 0 0];markerface=[1 1 1];



[pred60lms, pred60lm, pred60s]=plotFooofHarmonics(compiledData_HAF,3,lb,ub);
[LMS, LM, S, LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_HAF,lb,ub, pred60lms, pred60lm, pred60s);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 1 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(2)
subplot(3,2,2)
hold on
errorbar(TemporalFrequency,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Sum of harmonics'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

subplot(3,2,4)
hold on
errorbar(TemporalFrequency,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

subplot(3,2,6)
hold on
errorbar(TemporalFrequency([1:3 5]),Sm([1:3 5]),Sm([1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm([1:3 5]),'ob','LineWidth',2)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[0 0.04];

VEP_harmonic=cat(3,LMS,LM,S);
for i=1:size(compiledData_HAF,1)
    temp=nanmedian(compiledData_HAF(i).vds,3);
    temp2=squeeze(VEP_harmonic(i,:,:))';
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            VDS_all=cat(2,VDS_all,temp(j,k));
            VEP_all=cat(2,VEP_all,temp2(j,k));
            group=cat(2,group,{'HAF'});
            TF=cat(2,TF,TemporalFrequency(k));
        
            switch j
                case 1
                    PRP=cat(2,PRP,{'LMS'});
                case 2
                    PRP=cat(2,PRP,{'LM'});
                case 3
                    PRP=cat(2,PRP,{'S'});
            end
        end
    end
end

VDS_anova=anovan(VDS_all,{group,TF,PRP},'model','interaction','varnames',{'group','Temporal frequency','channel'});
VEP_anova=anovan(VEP_all,{group,TF,PRP},'model','interaction','varnames',{'group','Temporal frequency','channel'});

figure(6)
plotWithErrorbars(-1,LMSm(end),LMSci(:,end),markerline,markeredge,markerface)
title(['Luminance 30 Hz, R squared=' num2str(Rmva_vepHAf(1,2)^2) ', p=' num2str(Pmva_vepHAf(1,2))])
ylabel('amplitude at stimulus frequency (mV)')
xlabel('Number of headache days in past 3 months')
ax=gca; ax.TickDir='out'; ax.Box='off'; ax.YLim=[-0.001 0.04]; ax.XLim=[-2 31];



%% local functions

function [pred60lms, pred60lm, pred60s]=plotFooofHarmonics(compiledData,fig_num,lb,ub)

lin_fun=@(x,xdata)x(1)+(x(2)*xdata);
counter_lms=0;
counter_lm=0;
counter_s=0;

    for x=1:5     
        for y=1:size(compiledData,1)
            fooof_peak_harmonics=compiledData(y).fooof_peak_harmonics;
            x_harmonics=cell2mat(compiledData(y).fooof_peak_harmonics_freq(1,x));
            LMS(y,:)=cell2mat(fooof_peak_harmonics(1,x));
            LM(y,:)=cell2mat(fooof_peak_harmonics(2,x));
            S(y,:)=cell2mat(fooof_peak_harmonics(3,x));
        end

        Bootstat=bootstrp(1000,@nanmedian,LMS,1);
        Bootstat=sort(Bootstat,1);
        LMSci=Bootstat([lb ub],:);
        LMSm=Bootstat(500,:);
        
        Bootstat=bootstrp(1000,@nanmedian,LM,1);
        Bootstat=sort(Bootstat,1);
        LMci=Bootstat([lb ub],:);
        LMm=Bootstat(500,:);
        
        Bootstat=bootstrp(1000,@nanmedian,S,1);
        Bootstat=sort(Bootstat,1);
        Sci=Bootstat([lb ub],:);
        Sm=Bootstat(500,:);
        
        figure(fig_num)
        subplot(5,1,x)
        hold on
        errorbar(x_harmonics,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'-ok','LineWidth',2)
        errorbar(x_harmonics,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'-or','LineWidth',2)
        errorbar(x_harmonics,Sm,Sm-Sci(1,:),Sci(2,:)-Sm,'-ob','LineWidth',2)
        plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
        title(['harmonics'])
        ylabel('amplitude at stimulus frequency (mV)')
        xlabel('Stimulus frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 95];
        ax.YLim=[-0.001 0.01];
        
        if x>3
            freq_fine=round(x_harmonics(1)):1:round(x_harmonics(end));
            Hz60=find(freq_fine==60);
            x0=[0.01 -0.001];

            counter_lms=counter_lms+1;
            counter_lm=counter_lm+1;
            counter_s=counter_s+1;

            fit_harm=lsqcurvefit(lin_fun,x0,x_harmonics,LMSm);
            temp=lin_fun(fit_harm,freq_fine);
            pred60lms(counter_lms)=temp(Hz60);
        
            fit_harm=lsqcurvefit(lin_fun,x0,x_harmonics,LMm);
            temp=lin_fun(fit_harm,freq_fine);
            pred60lm(counter_lm)=temp(Hz60);
            
            fit_harm=lsqcurvefit(lin_fun,x0,x_harmonics,Sm);
            temp=lin_fun(fit_harm,freq_fine);
            pred60s(counter_s)=temp(Hz60);
        end
        
        clear LMS LM S LMSm LMm Sm LMSci LMci Sci
        
    end
    
end


function [LMS, LM, S, LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData,lb,ub, pred60lms, pred60lm, pred60s)

    for x=1:5     
        for y=1:size(compiledData,1)
            fooof_peak_harmonics=compiledData(y).fooof_peak_harmonics;
            x_harmonics=cell2mat(compiledData(y).fooof_peak_harmonics_freq(1,x));

            LMS(y,x)=sum(cell2mat(fooof_peak_harmonics(1,x)));
            LM(y,x)=sum(cell2mat(fooof_peak_harmonics(2,x)));
            S(y,x)=sum(cell2mat(fooof_peak_harmonics(3,x)));
            if x==4
                LMS(y,x)=LMS(y,x)+pred60lms(1);
                LM(y,x)=LMS(y,x)+pred60lm(1);
                S(y,x)=LMS(y,x)+pred60s(1);
            end
            
            if x==5
                LMS(y,x)=LMS(y,x)+pred60lms(2);
                LM(y,x)=LM(y,x)+pred60lm(2);
                S(y,x)=S(y,x)+pred60s(2);
            end
        end
    end
        Bootstat=bootstrp(1000,@nanmedian,LMS,1);
        Bootstat=sort(Bootstat,1);
        LMSci=Bootstat([lb ub],:);
        LMSm=Bootstat(500,:);
        
        Bootstat=bootstrp(1000,@nanmedian,LM,1);
        Bootstat=sort(Bootstat,1);
        LMci=Bootstat([lb ub],:);
        LMm=Bootstat(500,:);

        Bootstat=bootstrp(1000,@nanmedian,S,1);
        Bootstat=sort(Bootstat,1);
        Sci=Bootstat([lb ub],:);
        Sm=Bootstat(500,:);

end

function [VDS]=calcVDS(compiledData,lb,ub)
    
    LMS=[];
    LM=[];
    S=[];

    for x=1:size(compiledData,1)
        temp=compiledData(x).vds;
        LMS=cat(1,LMS,temp(1,:,:));
        LM=cat(1,LM,temp(2,:,:));
        S=cat(1,S,temp(3,:,:));
    end
    
    LMS=squeeze(nanmedian(LMS,3));
    Bootstat=bootstrp(1000,@nanmedian,LMS,1);
    Bootstat=sort(Bootstat,1);
    LMSCI=Bootstat([lb ub],:);
    LMSM=Bootstat(500,:);
    
    LM=squeeze(nanmedian(LM,3));
    Bootstat=bootstrp(1000,@nanmedian,LM,1);
    Bootstat=sort(Bootstat,1);
    LMCI=Bootstat([lb ub],:);
    LMM=Bootstat(500,:);
    
    S=squeeze(nanmedian(S,3));
    Bootstat=bootstrp(1000,@nanmedian,S,1);
    Bootstat=sort(Bootstat,1);
    SCI=Bootstat([lb ub],:);
    SM=Bootstat(500,:);
   
    VDS.LMSM=LMSM;
    VDS.LMSCI=LMSCI;
    VDS.LMM=LMM;
    VDS.LMCI=LMCI;
    VDS.SM=SM;
    VDS.SCI=SCI;

end

function [ttf_fit, TemporalFrequency_fit]=getTTFfits(VEPresponse,stimulusFreqHz,x0)
    % Fit the Watson model to data
    % Adjust the VEP response to deal with negative values
    minVEP=min(VEPresponse);
    if minVEP<0
        scaledVEP=VEPresponse-minVEP;
    else
        scaledVEP=VEPresponse;
        minVEP=0;
    end
    % Find the maximum interpolated VEP response
    stimulusFreqHzFine=logspace(0,log10(max(stimulusFreqHz)),100);
    splineInterpolatedMax=max(spline(stimulusFreqHz,scaledVEP,stimulusFreqHzFine));
    
    % Scale the x vector so that the max is 0
    scaledVEP=scaledVEP./splineInterpolatedMax;
    x0=cat(2,x0,max(scaledVEP));
    myObj=@(p)sqrt(sum((scaledVEP-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    x0_lb=x0-[0.5 1 0.5 0.001];
    x0_ub=x0+[0.5 1 0.5 0.001];
    params=fmincon(myObj,x0,[],[],[],[],x0_lb,x0_ub);
    
    
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP;
    TemporalFrequency_fit=stimulusFreqHzFine;
end

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
    'MELA_0209';'MELA_0204';...
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

        if compiledData(1).group=='MWVA'
            counter_MVA=counter_MVA+1;
            compiledData_MVA(counter_MVA,:)=compiledData;
        else if compiledData(1).group=='HA f'
                counter_HAF=counter_HAF+1;
               compiledData_HAF(counter_HAF,:)=compiledData;
            else
                disp('error: group');
            end
        end

end

clear ans compiledData counter_HAF counter_MVA observerID x filenameComp

% flicker discomfort
VDS=calcVDS(compiledData_MVA,compiledData_HAF,lb,ub);

% % Plot median psd for stimulus frequency across groups
% [LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_all,lb,ub);
% 
% [ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
% [ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[2 2 1]);
% [ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);
% 
% figure(1)
% subplot(1,4,1)
% hold on
% fillcolor=[0.95 0.95 0.95];edgecolor=[0.85 0.85 0.85];markeredge=[0.5 0.5 0.5];markerface=[0.5 0.5 0.5];
% plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
% plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
% ylabel('amplitude at stimulus frequency (mV)')
% title(['all subjects (n=' num2str(size(compiledData_all,1)) '), Background'])
% ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.001 0.02];
% 
% subplot(1,4,2)
% hold on
% fillcolor=[0.85 0.85 0.85];edgecolor=[0.75 0.75 0.75];markeredge=[0 0 0];markerface=[0 0 0];
% plotWithErrorfill(TemporalFrequency,LMSm,LMSci,edgecolor,fillcolor,markeredge,markerface)
% plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
% plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
% title(['LMS'])
% ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.001 0.02];
% 
% subplot(1,4,3)
% hold on
% fillcolor=[1 0.9 0.9];edgecolor=[1 0.8 0.8];markeredge=[1 0 0];markerface=[1 0 0];
% plotWithErrorfill(TemporalFrequency,LMm,LMci,edgecolor,fillcolor,markeredge,markerface)
% plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
% plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
% xlabel('Stimulus frequency')
% title(['LM'])
% ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.001 0.02];
% 
% 
% subplot(1,4,4)
% hold on
% fillcolor=[0.9 0.9 1];edgecolor=[0.8 0.8 1];markeredge=[0 0 1];markerface=[0 0 1];
% plotWithErrorfill(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),edgecolor,fillcolor,markeredge,markerface)
% plot(TemporalFrequency_fitS,ttf_fitS,'-b')
% plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
% title(['S'])
% ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.001 0.02];
% 
% 
% clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci

% % Plot the harmonics across subjects
% [pred60lms, pred60lm, pred60s]=plotFooofHarmonics(compiledData_all,3,lb,ub);
% [pred60lmsMVA, pred60lmMVA, pred60sMVA]=plotFooofHarmonics(compiledData_MVA,4,lb,ub);
% [pred60lmsHAF, pred60lmHAF, pred60sHAF]=plotFooofHarmonics(compiledData_HAF,5,lb,ub);

% % Plot the sum of harmonics across subjects
% sumFooofHarmonics(compiledData_all,10,TemporalFrequency,lb,ub, pred60lms, pred60lm, pred60s);
% [LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_MVA,11,TemporalFrequency,lb,ub, pred60lmsMVA, pred60lmMVA, pred60sMVA);
% [LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_HAF,12,TemporalFrequency,lb,ub, pred60lms, pred60lm, pred60s);


% Plot median psd for stimulus frequency between groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_MVA,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[2 2 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(2)
subplot(2,4,1)
hold on
fillcolor=[0.95 0.95 0.95];edgecolor=[0.85 0.85 0.85];markeredge=[0.5 0.5 0.5];markerface=[0.5 0.5 0.5];
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
title(['Migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) '), Background'])
ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.001 0.02];


subplot(2,4,2)
hold on
fillcolor=[0.85 0.85 0.85];edgecolor=[0.75 0.75 0.75];markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorfill(TemporalFrequency,LMSm,LMSci,edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];


subplot(2,4,3)
hold on
fillcolor=[1 0.9 0.9];edgecolor=[1 0.8 0.8];markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorfill(TemporalFrequency,LMm,LMci,edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];


subplot(2,4,4)
hold on
fillcolor=[0.9 0.9 1];edgecolor=[0.8 0.8 1];markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorfill(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];


% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
figure(7)
subplot(1,2,1)
hold on
markeredge=[0 0 0];markerface=[0 0 0];
plotWithXYErrorbars(LMSm,VDS.LMS_mvaM,LMSci,VDS.LMS_mvaCI,markeredge,markerface)
markeredge=[1 0 0];markerface=[1 0 0];
plotWithXYErrorbars(LMm,VDS.LM_mvaM,LMci,VDS.LM_mvaCI,markeredge,markerface)
markeredge=[0 0 1];markerface=[0 0 1];
plotWithXYErrorbars(Sm(:,[1:3 5]),VDS.S_mvaM(:,[1:3 5]),Sci(:,[1:3 5]),VDS.S_mvaCI(:,[1:3 5]),markeredge,markerface)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
[Rmva_vdsvep,Pmva_vdsvep]=corrcoef(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMS_mvaM,VDS.LM_mvaM,VDS.S_mvaM([1:3 5])));
title(['n=14, R squared=' num2str(Rmva_vdsvep(1,2)^2) ' , p=' num2str(Pmva_vdsvep(1,2))])
ylabel('flicker discomfort')
xlabel('amplitude at stimulus frequency (mV)')
ax=gca; ax.XLim=[-0.001 0.022]; ax.YLim=[0 10];


clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci



[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_HAF,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[2 2 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);
    
figure(2)
subplot(2,4,5)
hold on
fillcolor=[0.95 0.95 0.95];edgecolor=[0.85 0.85 0.85];markeredge=[0.5 0.5 0.5];markerface=[0.5 0.5 0.5];
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
title(['Headache free controls (n=' num2str(size(compiledData_HAF,1)) '), Background'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];

subplot(2,4,6)
hold on
fillcolor=[0.85 0.85 0.85];edgecolor=[0.75 0.75 0.75];markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorfill(TemporalFrequency,LMSm,LMSci,edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];

subplot(2,4,7)
hold on
fillcolor=[1 0.9 0.9];edgecolor=[1 0.8 0.8];markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorfill(TemporalFrequency,LMm,LMci,edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];

subplot(2,4,8)
hold on
fillcolor=[0.9 0.9 1];edgecolor=[0.8 0.8 1];markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorfill(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),edgecolor,fillcolor,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.001 0.02];

% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
figure(7)
subplot(1,2,2)
hold on
markeredge=[0 0 0];markerface=[1 1 1];
plotWithXYErrorbars(LMSm,VDS.LMS_hafM,LMSci,VDS.LMS_hafCI,markeredge,markerface)
markeredge=[1 0 0];markerface=[1 1 1];
plotWithXYErrorbars(LMm,VDS.LM_hafM,LMci,VDS.LM_hafCI,markeredge,markerface)
markeredge=[0 0 1];markerface=[1 1 1];
plotWithXYErrorbars(Sm(:,[1:3 5]),VDS.S_hafM(:,[1:3 5]),Sci(:,[1:3 5]),VDS.S_hafCI(:,[1:3 5]),markeredge,markerface)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
[Rhaf_vdsvep,Phaf_vdsvep]=corrcoef(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMS_hafM,VDS.LM_hafM,VDS.S_hafM([1:3 5])));
title(['n=14, R squared=' num2str(Rhaf_vdsvep(1,2)^2) ' , p=' num2str(Phaf_vdsvep(1,2))])
ylabel('flicker discomfort')
xlabel('amplitude at stimulus frequency (mV)')
ax=gca; ax.TickDir='out'; ax.Box='off'; ax.XLim=[-0.001 0.022]; ax.YLim=[0 10];

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci

% plot luminance 30 Hz VEP response as a function of headache frequency
figure(20)
for i=1:size(compiledData_MVA,1)
    X_mva(i)=table2array(compiledData_MVA(i).subject(:,4));
    Y_mva(i)=compiledData_MVA(i).fooof_peak_Fr(1,5);
end

for i=1:size(compiledData_HAF,1)
    X_haf(i)=table2array(compiledData_HAF(i).subject(:,4));
    Y_haf(i)=compiledData_HAF(i).fooof_peak_Fr(1,5);
end

[R_vepHAf,P_vepHAf]=corrcoef(cat(2,X_mva,Y_haf),cat(2,Y_mva,Y_haf));
[Rmva_vepHAf,Pmva_vepHAf]=corrcoef(X_mva,Y_mva);

plot(X_mva,Y_mva,'ok','MarkerFaceColor','k')
hold on
plot(-1*ones(size(X_haf)),Y_haf,'ok','MarkerFaceColor','w')
title(['Luminance 30 Hz, R squared=' num2str(Rmva_vepHAf(1,2)^2) ', p=' num2str(Pmva_vepHAf(1,2))])
ylabel('amplitude at stimulus frequency (mV)')
xlabel('Number of headache days in past 3 months')
ax=gca; ax.TickDir='out'; ax.Box='off'; ax.YLim=[-0.001 0.025]; ax.XLim=[-1 31];

% compare subject demographics across groups

%% local functions

function [LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData,lb,ub)
LMS=[];
LM=[];
S=[];

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
    
    Bootstat=bootstrp(1000,@nanmedian,BKGD,1);
    Bootstat=sort(Bootstat,1);
    BKGDci=Bootstat([lb ub],:);
    BKGDm=Bootstat(500,:);
end

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
        
        clear LMS LM S LMSm LMm Sm BKGDm LMSci LMci Sci
        
    end
    
end


function [LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData,fig_num,TemporalFrequency,lb,ub, pred60lms, pred60lm, pred60s)

    for x=1:5     
        for y=1:size(compiledData,1)
            fooof_peak_harmonics=compiledData(y).fooof_peak_harmonics;
            x_harmonics=cell2mat(compiledData(y).fooof_peak_harmonics_freq(1,x));
% % 1st-4th harmonic
%             LMS(y,x)=sum(cell2mat(fooof_peak_harmonics(1,x)));
%             LM(y,x)=sum(cell2mat(fooof_peak_harmonics(2,x)));
%             S(y,x)=sum(cell2mat(fooof_peak_harmonics(3,x)));
%             if x==4
%                 LMS(y,x)=LMS(y,x)+pred60lms(1);
%                 LM(y,x)=LMS(y,x)+pred60lm(1);
%                 S(y,x)=LMS(y,x)+pred60s(1);
%             end
%             
%             if x==5
%                 LMS(y,x)=LMS(y,x)+pred60lms(2);
%                 LM(y,x)=LM(y,x)+pred60lm(2);
%                 S(y,x)=S(y,x)+pred60s(2);
%             end

%1st and 2nd harmonic
            temp=cell2mat(fooof_peak_harmonics(1,x));
            LMS(y,x)=sum(temp(1:2));
            temp=cell2mat(fooof_peak_harmonics(2,x));
            LM(y,x)=sum(temp(1:2));
            temp=cell2mat(fooof_peak_harmonics(3,x));
            S(y,x)=sum(temp(1:2));
            
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
        
        figure(fig_num)
        hold on
        errorbar(TemporalFrequency,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'-ok','LineWidth',2)
        errorbar(TemporalFrequency,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'-or','LineWidth',2)
        errorbar(TemporalFrequency([1:3 5]),Sm([1:3 5]),Sm([1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm([1:3 5]),'-ob','LineWidth',2)
        plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
        title(['Sum of harmonics'])
        ylabel('amplitude at stimulus frequency (mV)')
        xlabel('Stimulus frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 95];
        ax.YLim=[-0.001 0.03];

end

function [VDS]=calcVDS(compiledData_MVA,compiledData_HAF,lb,ub)
    
    LMS_mva=[];
    LM_mva=[];
    S_mva=[];

    for x=1:size(compiledData_MVA,1)
        temp=compiledData_MVA(x).vds;
        LMS_mva=cat(1,LMS_mva,temp(1,:,:));
        LM_mva=cat(1,LM_mva,temp(2,:,:));
        S_mva=cat(1,S_mva,temp(3,:,:));
    end
    
    LMS_mva=squeeze(nanmedian(LMS_mva,3));
    Bootstat=bootstrp(1000,@nanmedian,LMS_mva,1);
    Bootstat=sort(Bootstat,1);
    LMS_mvaCI=Bootstat([lb ub],:);
    LMS_mvaM=Bootstat(500,:);
    
    LM_mva=squeeze(nanmedian(LM_mva,3));
    Bootstat=bootstrp(1000,@nanmedian,LM_mva,1);
    Bootstat=sort(Bootstat,1);
    LM_mvaCI=Bootstat([lb ub],:);
    LM_mvaM=Bootstat(500,:);
    
    S_mva=squeeze(nanmedian(S_mva,3));
    Bootstat=bootstrp(1000,@nanmedian,S_mva,1);
    Bootstat=sort(Bootstat,1);
    S_mvaCI=Bootstat([lb ub],:);
    S_mvaM=Bootstat(500,:);
    
    LMS_haf=[];
    LM_haf=[];
    S_haf=[];

    for x=1:size(compiledData_HAF,1)
        temp=compiledData_HAF(x).vds;
        LMS_haf=cat(1,LMS_haf,temp(1,:,:));
        LM_haf=cat(1,LM_haf,temp(2,:,:));
        S_haf=cat(1,S_haf,temp(3,:,:));
    end
    
    LMS_haf=squeeze(nanmedian(LMS_haf,3));
    Bootstat=bootstrp(1000,@nanmedian,LMS_haf,1);
    Bootstat=sort(Bootstat,1);
    LMS_hafCI=Bootstat([lb ub],:);
    LMS_hafM=Bootstat(500,:);
    
    LM_haf=squeeze(nanmedian(LM_haf,3));
    Bootstat=bootstrp(1000,@nanmedian,LM_haf,1);
    Bootstat=sort(Bootstat,1);
    LM_hafCI=Bootstat([lb ub],:);
    LM_hafM=Bootstat(500,:);
    
    S_haf=squeeze(nanmedian(S_haf,3));
    Bootstat=bootstrp(1000,@nanmedian,S_haf,1);
    Bootstat=sort(Bootstat,1);
    S_hafCI=Bootstat([lb ub],:);
    S_hafM=Bootstat(500,:);
    
    VDS.LMS_mvaM=LMS_mvaM;
    VDS.LMS_mvaCI=LMS_mvaCI;
    VDS.LM_mvaM=LM_mvaM;
    VDS.LM_mvaCI=LM_mvaCI;
    VDS.S_mvaM=S_mvaM;
    VDS.S_mvaCI=S_mvaCI;

    VDS.LMS_hafM=LMS_hafM;
    VDS.LMS_hafCI=LMS_hafCI;
    VDS.LM_hafM=LM_hafM;
    VDS.LM_hafCI=LM_hafCI;
    VDS.S_hafM=S_hafM;
    VDS.S_hafCI=S_hafCI;

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
%     params=lsqcurvefit(myObj,x0,stimulusFreqHz,scaledVEP);
    
    figure(100)
    semilogx(stimulusFreqHzFine,watsonTemporalModelvep(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP,'-k');
    hold on
    semilogx(stimulusFreqHz, VEPresponse, '*r');
    hold off
    
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP;
    TemporalFrequency_fit=stimulusFreqHzFine;
end

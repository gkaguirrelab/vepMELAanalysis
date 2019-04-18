% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169';...
    'MELA_0194';'MELA_0179';...
    'MELA_0191';'MELA_0174';...
    'MELA_0120'];

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
        else if compiledData(1).group=='HA f'
                counter_HAF=counter_HAF+1;
               compiledData_HAF(counter_HAF,:)=compiledData;
            else
                disp('error: group');
            end
        end

end

clear ans compiledData counter_HAF counter_MVA observerID x filenameComp

% visual discomfort scale
VDS=calcVDS(compiledData_MVA,compiledData_HAF,lb,ub);

% Plot median psd for stimulus frequency across groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_all,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[4 2 1 max(LMSm)]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 2 1 max(LMm)]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[4 2 1 max(Sm)]);

figure(1)
subplot(1,4,1)
hold on
X=cat(2,TemporalFrequency,fliplr(TemporalFrequency));
CI=cat(2,BKGDci(1,:),fliplr(BKGDci(2,:)));
fill(X,CI,[0.95 0.95 0.95],'EdgeColor',[0.85 0.85 0.85]);
plot(TemporalFrequency,BKGDm,'ok','MarkerFaceColor',[0.5 0.5 0.5])
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('power spectra at stimulus frequency')
title(['all subjects (n=' num2str(size(compiledData_all,1)) '), Background'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,2)
hold on
CI=cat(2,LMSci(1,:),fliplr(LMSci(2,:)));
fill(X,CI,[0.9 0.9 0.9],'EdgeColor',[0.75 0.75 0.75]);
plot(TemporalFrequency,LMSm,'ok','MarkerFaceColor','k')
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['LMS'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,3)
hold on
CI=cat(2,LMci(1,:),fliplr(LMci(2,:)));
fill(X,CI,[1 0.9 0.9],'EdgeColor',[1 0.8 0.8]);
plot(TemporalFrequency,LMm,'or','MarkerFaceColor','r')
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,4)
hold on
X=cat(2,TemporalFrequency([1:3 5]),fliplr(TemporalFrequency([1:3 5])));
CI=cat(2,Sci(1,[1:3 5]),fliplr(Sci(2,[1:3 5])));
fill(X,CI,[0.9 0.9 1],'EdgeColor',[0.8 0.8 1]);
plot(TemporalFrequency([1:3 5]),Sm([1:3 5]),'ob','MarkerFaceColor','b')
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['S'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci

% Plot the harmonics across subjects
[lmsFit_harmonics, lmFit_harmonics, Sfit_harmonics]=plotFooofHarmonics(compiledData_all,3,lb,ub);
[lmsFit_harmonicsMVA, lmFit_harmonicsMVA, Sfit_harmonicsMVA]=plotFooofHarmonics(compiledData_MVA,4,lb,ub);
[lmsFit_harmonicsHAF, lmFit_harmonicsHAF, Sfit_harmonicsHAF]=plotFooofHarmonics(compiledData_HAF,5,lb,ub);

% Plot the sum of harmonics across subjects
sumFooofHarmonics(compiledData_all,10,TemporalFrequency,lb,ub);
[LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_MVA,11,TemporalFrequency,lb,ub);
% plot visual discomfort data as a function of VEP power at the stimulus
% frequency+harmonics
figure(7)
subplot(2,2,3)
hold on
errorbar(LMSm,VDS.LMS_mvaM,VDS.LMS_mvaM-VDS.LMS_mvaCI(1,:),VDS.LMS_mvaCI(2,:)-VDS.LMS_mvaM,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2,'MarkerFaceColor','k')
errorbar(LMm,VDS.LM_mvaM,VDS.LM_mvaM-VDS.LM_mvaCI(1,:),VDS.LM_mvaCI(2,:)-VDS.LM_mvaM,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2,'MarkerFaceColor','r')
errorbar(Sm(:,[1:3 5]),VDS.S_mvaM(:,[1:3 5]),VDS.S_mvaM(:,[1:3 5])-VDS.S_mvaCI(1,[1:3 5]),VDS.S_mvaCI(2,[1:3 5])-VDS.S_mvaM(:,[1:3 5]),...
    Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) ') fundamental+harmonics'])
ylabel('visual discomfort scale')
xlabel('VEP power at stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XLim=[-0.001 0.021];
ax.YLim=[0 10];
mvaFit_vdsvep_harm=fitlm(cat(2,LMSm,LMm,Sm),cat(2,VDS.LMS_mvaM,VDS.LM_mvaM,VDS.S_mvaM));


[LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData_HAF,12,TemporalFrequency,lb,ub);
figure(7)
subplot(2,2,4)
hold on
errorbar(LMSm,VDS.LMS_hafM,VDS.LMS_hafM-VDS.LMS_hafCI(1,:),VDS.LMS_hafCI(2,:)-VDS.LMS_hafM,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2,'MarkerFaceColor','w')
errorbar(LMm,VDS.LM_hafM,VDS.LM_hafM-VDS.LM_hafCI(1,:),VDS.LM_hafCI(2,:)-VDS.LM_hafM,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2,'MarkerFaceColor','w')
errorbar(Sm(:,[1:3 5]),VDS.S_hafM(:,[1:3 5]),VDS.S_hafM(:,[1:3 5])-VDS.S_hafCI(1,[1:3 5]),VDS.S_hafCI(2,[1:3 5])-VDS.S_hafM(:,[1:3 5]),...
    Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','w')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Headache free control (n=' num2str(size(compiledData_HAF,1)) ') fundamental+harmonics'])
ylabel('visual discomfort scale')
xlabel('VEP power at stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XLim=[-0.001 0.021];
ax.YLim=[0 10];
hafFit_vdsvep_harm=fitlm(cat(2,LMSm,LMm,Sm),cat(2,VDS.LMS_hafM,VDS.LM_hafM,VDS.S_hafM));


% Plot median psd for stimulus frequency between groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_MVA,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[4 2 1 max(LMSm)]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 2 1 max(LMm)]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[4 2 1 max(Sm)]);

figure(2)
subplot(1,4,1)
hold on
X=cat(2,TemporalFrequency,fliplr(TemporalFrequency));
CI=cat(2,BKGDci(1,:),fliplr(BKGDci(2,:)));
fill(X,CI,[0.95 0.95 0.95],'EdgeColor',[0.85 0.85 0.85]);
plot(TemporalFrequency,BKGDm,'ok','MarkerFaceColor',[0.5 0.5 0.5])
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('power spectra at stimulus frequency')
title(['migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) '), Background'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,2)
hold on
CI=cat(2,LMSci(1,:),fliplr(LMSci(2,:)));
fill(X,CI,[0.9 0.9 0.9],'EdgeColor',[0.75 0.75 0.75]);
plot(TemporalFrequency,LMSm,'ok','MarkerFaceColor','k')
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['LMS'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,3)
hold on
CI=cat(2,LMci(1,:),fliplr(LMci(2,:)));
fill(X,CI,[1 0.9 0.9],'EdgeColor',[1 0.8 0.8]);
plot(TemporalFrequency,LMm,'or','MarkerFaceColor','r')
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

subplot(1,4,4)
hold on
X=cat(2,TemporalFrequency([1:3 5]),fliplr(TemporalFrequency([1:3 5])));
CI=cat(2,Sci(1,[1:3 5]),fliplr(Sci(2,[1:3 5])));
fill(X,CI,[0.9 0.9 1],'EdgeColor',[0.8 0.8 1]);
plot(TemporalFrequency([1:3 5]),Sm([1:3 5]),'ob','MarkerFaceColor','b')
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['S'])
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
figure(7)
subplot(2,2,1)
hold on
errorbar(LMSm,VDS.LMS_mvaM,VDS.LMS_mvaM-VDS.LMS_mvaCI(1,:),VDS.LMS_mvaCI(2,:)-VDS.LMS_mvaM,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2,'MarkerFaceColor','k')
errorbar(LMm,VDS.LM_mvaM,VDS.LM_mvaM-VDS.LM_mvaCI(1,:),VDS.LM_mvaCI(2,:)-VDS.LM_mvaM,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2,'MarkerFaceColor','r')
errorbar(Sm(:,[1:3 5]),VDS.S_mvaM(:,[1:3 5]),VDS.S_mvaM(:,[1:3 5])-VDS.S_mvaCI(1,[1:3 5]),VDS.S_mvaCI(2,[1:3 5])-VDS.S_mvaM(:,[1:3 5]),...
    Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) ')'])
ylabel('visual discomfort scale')
xlabel('VEP power at stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XLim=[-0.001 0.021];
ax.YLim=[0 10];

mvaFit_vdsvep=fitlm(cat(2,LMSm,LMm,Sm),cat(2,VDS.LMS_mvaM,VDS.LM_mvaM,VDS.S_mvaM));

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci



[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_HAF,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[4 2 1 max(LMSm)]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 2 1 max(LMm)]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[4 2 1 max(Sm)]);
    
figure(2)
subplot(1,4,1)
hold on
X=cat(2,TemporalFrequency,fliplr(TemporalFrequency));
CI=cat(2,BKGDci(1,:),fliplr(BKGDci(2,:)));
fill(X,CI,[0.95 0.95 0.95],'EdgeColor',[0.85 0.85 0.85]);
plot(TemporalFrequency,BKGDm,'ok','MarkerFaceColor',[0.5 0.5 0.5])
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])

subplot(1,4,2)
hold on
CI=cat(2,LMSci(1,:),fliplr(LMSci(2,:)));
fill(X,CI,[0.9 0.9 0.9],'EdgeColor',[0.75 0.75 0.75]);
plot(TemporalFrequency,LMSm,'ok','MarkerFaceColor','k')
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['headache free (n=' num2str(size(compiledData_HAF,1)) '), LMS'])


subplot(1,4,3)
hold on
CI=cat(2,LMci(1,:),fliplr(LMci(2,:)));
fill(X,CI,[1 0.9 0.9],'EdgeColor',[1 0.8 0.8]);
plot(TemporalFrequency,LMm,'or','MarkerFaceColor','r')
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])


subplot(1,4,4)
hold on
X=cat(2,TemporalFrequency([1:3 5]),fliplr(TemporalFrequency([1:3 5])));
CI=cat(2,Sci(1,[1:3 5]),fliplr(Sci(2,[1:3 5])));
fill(X,CI,[0.9 0.9 1],'EdgeColor',[0.8 0.8 1]);
plot(TemporalFrequency([1:3 5]),Sm([1:3 5]),'ob','MarkerFaceColor','b')
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])

% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
figure(7)
subplot(2,2,2)
hold on
errorbar(LMSm,VDS.LMS_hafM,VDS.LMS_hafM-VDS.LMS_hafCI(1,:),VDS.LMS_hafCI(2,:)-VDS.LMS_hafM,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2,'MarkerFaceColor','w')
errorbar(LMm,VDS.LM_hafM,VDS.LM_hafM-VDS.LM_hafCI(1,:),VDS.LM_hafCI(2,:)-VDS.LM_hafM,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2,'MarkerFaceColor','w')
errorbar(Sm(:,[1:3 5]),VDS.S_hafM(:,[1:3 5]),VDS.S_hafM(:,[1:3 5])-VDS.S_hafCI(1,[1:3 5]),VDS.S_hafCI(2,[1:3 5])-VDS.S_hafM(:,[1:3 5]),...
    Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','w')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Headache free control (n=' num2str(size(compiledData_HAF,1)) ')'])
ylabel('visual discomfort scale')
xlabel('VEP power at stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XLim=[-0.001 0.021];
ax.YLim=[0 10];

hafFit_vdsvep=fitlm(cat(2,LMSm,LMm,Sm),cat(2,VDS.LMS_hafM,VDS.LM_hafM,VDS.S_hafM));


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

Fit_vepHAfr=fitlm(cat(2,X_mva,Y_haf),cat(2,Y_mva,Y_haf));
mvaFit_vepHAfr=fitlm(X_mva,Y_mva);

plot(X_mva,Y_mva,'ok','MarkerFaceColor','k')
hold on
plot(X_haf,Y_haf,'ok','MarkerFaceColor','w')
title(['Luminance 30 Hz'])
ylabel('power at stimulus frequency')
xlabel('Number of headache days in past 3 months')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.YLim=[-0.001 0.012];
ax.XLim=[-1 31];

%% local functions

function [LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData,lb,ub)
LMS=[];
LM=[];
S=[];
BKGD=[];

    for x=1:size(compiledData,1)
        temp=compiledData(x).fooof_peak_Fr;
        temp2=compiledData(x).fooof_bkgd_Fr;
        LMS=cat(1,LMS,temp(1,:));
        LM=cat(1,LM,temp(2,:));
        S=cat(1,S,temp(3,:));
        BKGD=cat(1,BKGD,temp2);
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
    
    Bootstat=bootstrp(1000,@nanmedian,BKGD,1);
    Bootstat=sort(Bootstat,1);
    BKGDci=Bootstat([lb ub],:);
    BKGDm=Bootstat(500,:);
end

function [lmsFit_harmonics, lmFit_harmonics, sFit_harmonics]=plotFooofHarmonics(compiledData,fig_num,lb,ub)

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
        lmsFit_harmonics{x,:}=fitlm(x_harmonics,LMSm);
        
        Bootstat=bootstrp(1000,@nanmedian,LM,1);
        Bootstat=sort(Bootstat,1);
        LMci=Bootstat([lb ub],:);
        LMm=Bootstat(500,:);
        lmFit_harmonics{x,:}=fitlm(x_harmonics,LMm);

        Bootstat=bootstrp(1000,@nanmedian,S,1);
        Bootstat=sort(Bootstat,1);
        Sci=Bootstat([lb ub],:);
        Sm=Bootstat(500,:);
        sFit_harmonics{x,:}=fitlm(x_harmonics,Sm);
        
        figure(fig_num)
        subplot(5,1,x)
        hold on
        errorbar(x_harmonics,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'-ok','LineWidth',2)
        errorbar(x_harmonics,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'-or','LineWidth',2)
        errorbar(x_harmonics,Sm,Sm-Sci(1,:),Sci(2,:)-Sm,'-ob','LineWidth',2)
        plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
        title(['harmonics'])
        ylabel('power spectra at stimulus frequency')
        xlabel('Stimulus frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 95];
        ax.YLim=[-0.001 0.01];

        clear LMS LM S LMSm LMm Sm BKGDm LMSci LMci Sci
    end
end


function [LMSm, LMm, Sm, LMSci, LMci, Sci]=sumFooofHarmonics(compiledData,fig_num,TemporalFrequency,lb,ub)

    for x=1:5     
        for y=1:size(compiledData,1)
            fooof_peak_harmonics=compiledData(y).fooof_peak_harmonics;
            x_harmonics=cell2mat(compiledData(y).fooof_peak_harmonics_freq(1,x));
            LMS(y,x)=sum(cell2mat(fooof_peak_harmonics(1,x)));
            LM(y,x)=sum(cell2mat(fooof_peak_harmonics(2,x)));
            S(y,x)=sum(cell2mat(fooof_peak_harmonics(3,x)));
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
        errorbar(TemporalFrequency,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2)
        errorbar(TemporalFrequency,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2)
        errorbar(TemporalFrequency([1:3 5]),Sm([1:3 5]),Sm([1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm([1:3 5]),'ob','LineWidth',2)
        plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
        title(['Sum of harmonics'])
        ylabel('power spectra at stimulus frequency')
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
    stimulusFreqHzFine = logspace(0,log10(max(stimulusFreqHz)),100);
    splineInterpolatedMax = max(spline(stimulusFreqHz,scaledVEP,stimulusFreqHzFine));
    % Scale the x vector so that the max is zero
    scaledVEP=scaledVEP./splineInterpolatedMax;
    myObj=@(p)sqrt(sum((scaledVEP-watsonTemporalModel(stimulusFreqHz,p)).^2));
    params = fmincon(myObj,x0,[],[]);
    figure(100)
    semilogx(stimulusFreqHzFine,watsonTemporalModel(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP,'-k');
    hold on
    semilogx(stimulusFreqHz, VEPresponse, '*r');
    hold off
    
    ttf_fit=watsonTemporalModel(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP;
    TemporalFrequency_fit=stimulusFreqHzFine;
end

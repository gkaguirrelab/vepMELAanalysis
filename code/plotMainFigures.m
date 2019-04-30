% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169';...
    'MELA_0194';'MELA_0179';...
    'MELA_0191';'MELA_0174';...
    'MELA_0120';'MELA_0171';...
    ];

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

figure(1)
subplot(3,2,1)
hold on
markerline='-ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LMS_mvaM,VDS.LMS_mvaCI,markerline,markeredge,markerface)
title(['Flicker discomfort MwA'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,3)
hold on
markerline='-or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LM_mvaM,VDS.LM_mvaCI,markerline,markeredge,markerface)
ylabel(['FLicker dicomfort'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,5)
hold on
markerline='-ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),VDS.S_mvaM(1,[1:3 5]),VDS.S_mvaCI(:,[1:3 5]),markerline,markeredge,markerface)
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,2)
hold on
markerline='-ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LMS_hafM,VDS.LMS_hafCI,markerline,markeredge,markerface)
title(['Flicker discomfort HAf'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,4)
hold on
markerline='-or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LM_hafM,VDS.LM_hafCI,markerline,markeredge,markerface)
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,6)
hold on
markerline='-ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),VDS.S_hafM(1,[1:3 5]),VDS.S_hafCI(:,[1:3 5]),markerline,markeredge,markerface)
xlabel(['Stimulus frequency (Hz)'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];


% Plot median psd for stimulus frequency between groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_MVA,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[2 2 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(2)

markerline='ok';markeredge=[0.5 0.5 0.5];markerface=[0.5 0.5 0.5];
plotWithErrorbars(TemporalFrequency,BKGDm,BKGDci,markerline,markeredge,markerface)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
ylabel('amplitude at stimulus frequency (mV)')
title(['Migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) '), Background'])
ax=gca;ax.XScale='log'; ax.XLim=[0.95 35]; ax.YLim=[-0.002 0.022];


subplot(3,2,1)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none'; 
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,LMSm,LMSci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,2,3)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none';
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,LMm,LMci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,2,5)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none';
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


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

subplot(3,2,2)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none';
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='ok';markeredge=[0 0 0];markerface=[1 1 1];
plotWithErrorbars(TemporalFrequency,LMSm,LMSci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];

subplot(3,2,4)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none';
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='or';markeredge=[1 0 0];markerface=[1 1 1];
plotWithErrorbars(TemporalFrequency,LMm,LMci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];

subplot(3,2,6)
hold on
edgecolor='none';fillcolor=[0.8 0.8 0.8];markeredge='none';markerface='none';
plotWithErrorfill(TemporalFrequency,BKGDm,BKGDci,edgecolor,fillcolor,markeredge,markerface)
markerline='ob';markeredge=[0 0 1];markerface=[1 1 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];

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
    X_mva(i)=table2array(compiledData_MVA(i).subject(:,6));
    Y_mva(i)=compiledData_MVA(i).fooof_peak_Fr(1,5);
end

for i=1:size(compiledData_HAF,1)
    X_haf(i)=table2array(compiledData_HAF(i).subject(:,6));
    Y_haf(i)=compiledData_HAF(i).fooof_peak_Fr(1,5);
end

X_haf(find(isnan(X_haf)))=1;

[R_vepHAf,P_vepHAf]=corrcoef(cat(2,X_mva,Y_haf),cat(2,Y_mva,Y_haf));
[Rmva_vepHAf,Pmva_vepHAf]=corrcoef(X_mva,Y_mva);

plot(X_mva,Y_mva,'ok','MarkerFaceColor','k')
hold on
plot(-1*ones(size(X_haf)),Y_haf,'ok','MarkerFaceColor','w')
title(['Luminance 30 Hz, R squared=' num2str(Rmva_vepHAf(1,2)^2) ', p=' num2str(Pmva_vepHAf(1,2))])
ylabel('amplitude at stimulus frequency (mV)')
xlabel('Number of headache days in past 3 months')
ax=gca; ax.TickDir='out'; ax.Box='off'; ax.YLim=[-0.001 0.025]; ax.XLim=[-2 31];

% compare subject demographics across groups

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
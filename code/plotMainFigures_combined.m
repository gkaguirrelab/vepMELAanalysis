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
    'MELA_0209';'MELA_0204'];

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
        compiledData_ALL(x,:)=compiledData;

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
VDS=calcVDS(compiledData_ALL,lb,ub);
VDS_mva=calcVDS(compiledData_MVA,lb,ub);
VDS_haf=calcVDS(compiledData_HAF,lb,ub);

% calculate ANOVA for flicker discomfort by group, post-receptoral pathway,
% and temporal frequency

VEP_all=[];
VDS_all=[];
group=[];
PRP=[];
TF=[];
for i=1:size(compiledData_MVA,1)
    temp=nanmedian(compiledData_MVA(i).vds,3);
    temp2=compiledData_MVA(i).fooof_peak_Fr;
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

for i=1:size(compiledData_HAF,1)
    temp=nanmedian(compiledData_HAF(i).vds,3);
    temp2=compiledData_HAF(i).fooof_peak_Fr;
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

figure(100)
temp=1:120;
temp2=121:length(VDS_all);
errorbar(1,mean(VDS_all(temp)),std(VDS_all(temp))./sqrt(length(temp)-1),'ok')
hold on
errorbar(1,mean(VDS_all(temp2)),std(VDS_all(temp2))./sqrt(length(temp2)-1),'ob')

VDS_anova=anovan(VDS_all,{group,TF,PRP},'model','interaction','varnames',{'group','Temporal frequency','channel'});
VEP_anova=anovan(VEP_all,{group,TF,PRP},'model','interaction','varnames',{'group','Temporal frequency','channel'});

[ttf_fitLMS_vds,TemporalFrequency_fitLMS_vds]=getTTFfits(VDS.LMS_M,TemporalFrequency,[1 2 1]);
[ttf_fitLM_vds,TemporalFrequency_fitLM_vds]=getTTFfits(VDS.LM_M,TemporalFrequency,[2 2 1]);
[ttf_fitS_vds,TemporalFrequency_fitS_vds]=getTTFfits(VDS.S_M(1,[1:3 5]),TemporalFrequency([1:3 5]),[5.5 1 1]);

figure(3)
subplot(3,1,1)
hold on
markerline='ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LMS_M,VDS.LMS_CI,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS_vds,ttf_fitLMS_vds,'-k')
title(['Flicker discomfort MwA'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,1,2)
hold on
markerline='or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LM_M,VDS.LM_CI,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM_vds,ttf_fitLM_vds,'-r')
ylabel(['FLicker dicomfort'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,1,3)
hold on
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),VDS.S_M(1,[1:3 5]),VDS.S_CI(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS_vds,ttf_fitS_vds,'-b')
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];


% Plot median psd for stimulus frequency
[LMSm, LMm, Sm, LMSci, LMci, Sci]=medianFooofFrequency(compiledData_ALL,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM]=getTTFfits(LMm,TemporalFrequency,[4 1 1]);
[ttf_fitS,TemporalFrequency_fitS]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

figure(4)
subplot(3,1,1)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,LMSm,LMSci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,1,2)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,LMm,LMci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,1,3)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


% plot visual discomfort data as a function of VEP power at the stimulus
% frequency
figure(5)
hold on
plot(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMS_M,VDS.LM_M,VDS.S_M([1:3 5])),'.','Color',[0.5 0.5 0.5])
lsline
plot(LMSm,VDS.LMS_M,'ok','MarkerFaceColor','k','MarkerSize',8)
plot(LMm,VDS.LM_M,'or','MarkerFaceColor','r','MarkerSize',8)
plot(Sm([1:3 5]),VDS.S_M([1:3 5]),'ob','MarkerFaceColor','b','MarkerSize',8)
[R_vdsvep,P_vdsvep]=corrcoef(cat(2,LMSm,LMm,Sm([1:3 5])),cat(2,VDS.LMS_M,VDS.LM_M,VDS.S_M([1:3 5])));
title(['n=14, R squared=' num2str(R_vdsvep(1,2)^2) ' , p=' num2str(P_vdsvep(1,2))])
ylabel('flicker discomfort')
xlabel('amplitude at stimulus frequency (mV)')
ax=gca; ax.XLim=[-0.002 0.015]; ax.YLim=[0 10];


clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci


% plot luminance 30 Hz VEP response as a function of headache frequency
figure(6)
for i=1:size(compiledData_MVA,1)
    X_mva(i)=table2array(compiledData_MVA(i).subject(:,6));
    Y_mva(i)=compiledData_MVA(i).fooof_peak_Fr(1,5);
end

[Rmva_vepHAf,Pmva_vepHAf]=corrcoef(X_mva,Y_mva);
[LMSm, LMm, Sm, LMSci, LMci, Sci]=medianFooofFrequency(compiledData_HAF,lb,ub);
plot(X_mva,Y_mva,'ok','MarkerFaceColor','k','MarkerSize',8)
hold on
lsline
markerline='ok';markeredge=[0 0 0];markerface=[1 1 1];
plotWithErrorbars(-1,LMSm(end),LMSci(:,end),markerline,markeredge,markerface)
title(['Luminance 30 Hz, R squared=' num2str(Rmva_vepHAf(1,2)^2) ', p=' num2str(Pmva_vepHAf(1,2))])
ylabel('amplitude at stimulus frequency (mV)')
xlabel('Number of headache days in past 3 months')
ax=gca; ax.TickDir='out'; ax.Box='off'; ax.YLim=[-0.001 0.025]; ax.XLim=[-2 31];



%% local functions

function [LMSm, LMm, Sm, LMSci, LMci, Sci]=medianFooofFrequency(compiledData,lb,ub)
LMS=[];
LM=[];
S=[];

    for x=1:size(compiledData,1)
        temp=compiledData(x).fooof_peak_Fr;
        LMS=cat(1,LMS,temp(1,:));
        LM=cat(1,LM,temp(2,:));
        S=cat(1,S,temp(3,:));
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
    LMS_CI=Bootstat([lb ub],:);
    LMS_M=Bootstat(500,:);
    
    LM=squeeze(nanmedian(LM,3));
    Bootstat=bootstrp(1000,@nanmedian,LM,1);
    Bootstat=sort(Bootstat,1);
    LM_CI=Bootstat([lb ub],:);
    LM_M=Bootstat(500,:);
    
    S=squeeze(nanmedian(S,3));
    Bootstat=bootstrp(1000,@nanmedian,S,1);
    Bootstat=sort(Bootstat,1);
    S_CI=Bootstat([lb ub],:);
    S_M=Bootstat(500,:);
    
    VDS.LMS_M=LMS_M;
    VDS.LMS_CI=LMS_CI;
    VDS.LM_M=LM_M;
    VDS.LM_CI=LM_CI;
    VDS.S_M=S_M;
    VDS.S_CI=S_CI;

end
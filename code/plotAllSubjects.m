% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169';...
    'MELA_0194'];

counter_MVA=0;
counter_HAF=0;
A=[1.625 3.25 7.5 15 30];
ub=750;
lb=250;

savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

for x=1:size(subjects,1)
    
        observerID=subjects(x,:);
        savePath=fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
        'vepMELAanalysis','allChannels');

        filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);

        open(filenameComp);
        compiledData=ans.compiledData;
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


% Plot median psd for stimulus frequency across groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_all,lb,ub);

figure(1)
hold on
errorbar(A,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'k','LineWidth',2)
errorbar(A,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'r','LineWidth',2)
errorbar(A([1:3 5]),Sm(:,[1:3 5]),Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'b','LineWidth',2)
errorbar(A,BKGDm,BKGDm-BKGDci(1,:),BKGDci(2,:)-BKGDm,'Color',[0.5 0.5 0.5],'LineWidth',2)
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['all subjects (n=' num2str(size(compiledData_all,1)) ')'])
ylabel('power spectra at stimulus frequency')
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.01];

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci


% Plot median psd for stimulus frequency between groups
[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_MVA,lb,ub);

figure(2)
subplot(2,1,1)
hold on
errorbar(A,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'-ok','LineWidth',2,'MarkerFaceColor','k')
errorbar(A,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'-or','LineWidth',2,'MarkerFaceColor','r')
errorbar(A([1:3 5]),Sm(:,[1:3 5]),Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'-ob','LineWidth',2,'MarkerFaceColor','b')
errorbar(A,BKGDm,BKGDm-BKGDci(1,:),BKGDci(2,:)-BKGDm,'Color',[0.5 0.5 0.5],'LineWidth',2,'MarkerFaceColor',[0.5 0.5 0.5])
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Migraine with visual aura (n=' num2str(size(compiledData_MVA,1)) ')'])
ylabel('power spectra at stimulus frequency')
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.02];

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci



[LMSm, LMm, Sm, BKGDm, LMSci, LMci, Sci, BKGDci]=medianFooofFrequency(compiledData_HAF,lb,ub);

figure(2)
subplot(2,1,2)
hold on
errorbar(A,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'-ok','LineWidth',2,'MarkerFaceColor','w')
errorbar(A,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'-or','LineWidth',2,'MarkerFaceColor','w')
errorbar(A([1:3 5]),Sm(:,[1:3 5]),Sm(:,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(:,[1:3 5]),'-ob','LineWidth',2,'MarkerFaceColor','w')
errorbar(A,BKGDm,BKGDm-BKGDci(1,:),BKGDci(2,:)-BKGDm,'Color',[0.5 0.5 0.5],'LineWidth',2,'MarkerFaceColor','w')
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
title(['Headache free (n=' num2str(size(compiledData_HAF,1)) ')'])
ylabel('power spectra at stimulus frequency')
xlabel('Stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 35];
ax.YLim=[-0.001 0.02];

clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci

% Plot harmonics across all subjects
plotFooofHarmonics(compiledData_all,3,A,lb,ub);

plotFooofHarmonics(compiledData_MVA,4,A,lb,ub);

plotFooofHarmonics(compiledData_HAF,5,A,lb,ub);

% plot visual discomfort scale
plotVDS(compiledData_MVA,compiledData_HAF,6,A,lb,ub);

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



function []=plotFooofHarmonics(compiledData,fig_num,A,lb,ub)

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
        errorbar(x_harmonics,LMSm,LMSm-LMSci(1,:),LMSci(2,:)-LMSm,'ok','LineWidth',2)
        errorbar(x_harmonics,LMm,LMm-LMci(1,:),LMci(2,:)-LMm,'or','LineWidth',2)
        errorbar(x_harmonics,Sm,Sm-Sci(1,:),Sci(2,:)-Sm,'ob','LineWidth',2)
        plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
        title(['Stimulus frequency=' num2str(A(x))])
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

function []=plotVDS(compiledData_MVA,compiledData_HAF,fig_num,A,lb,ub)
    
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

    figure(fig_num)
   
    subplot(1,2,2)
    hold on
    errorbar(A,LMS_mvaM,LMS_mvaM-LMS_mvaCI(1,:),LMS_mvaCI(2,:)-LMS_mvaM,'-ok','LineWidth',2,'MarkerFaceColor','k')
    errorbar(A,LM_mvaM,LM_mvaM-LM_mvaCI(1,:),LM_mvaCI(2,:)-LM_mvaM,'-or','LineWidth',2,'MarkerFaceColor','r')
    errorbar(A([1:3 5]),S_mvaM(:,[1:3 5]),S_mvaM(:,[1:3 5])-S_mvaCI(1,[1:3 5]),S_mvaCI(2,[1:3 5])-S_mvaM(:,[1:3 5]),'-ob','LineWidth',2,'MarkerFaceColor','b')
    ylabel('visual discomfort scale')
    xlabel('temporal frequency of stimulus')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[0 11];

    subplot(1,2,1)
    hold on
    errorbar(A,LMS_hafM,LMS_hafM-LMS_hafCI(1,:),LMS_hafCI(2,:)-LMS_hafM,'-ok','LineWidth',2,'MarkerFaceColor','w')
    errorbar(A,LM_hafM,LM_hafM-LM_hafCI(1,:),LM_hafCI(2,:)-LM_hafM,'-or','LineWidth',2,'MarkerFaceColor','w')
    errorbar(A([1:3 5]),S_hafM(:,[1:3 5]),S_hafM(:,[1:3 5])-S_hafCI(1,[1:3 5]),S_hafCI(2,[1:3 5])-S_hafM(:,[1:3 5]),'-ob','LineWidth',2,'MarkerFaceColor','w')
    ylabel('visual discomfort scale')
    xlabel('temporal frequency of stimulus')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XScale='log';
    ax.XLim=[0.95 35];
    ax.YLim=[0 11];
end

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
ub=50;
lb=950;

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
errorbar(A([1:3 5]),Sm(1,[1:3 5]),Sm(1,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(1,[1:3 5]),'b','LineWidth',2)
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
errorbar(A([1:3 5]),Sm(1,[1:3 5]),Sm(1,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(1,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','b')
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
errorbar(A([1:3 5]),Sm(1,[1:3 5]),Sm(1,[1:3 5])-Sci(1,[1:3 5]),Sci(2,[1:3 5])-Sm(1,[1:3 5]),'ob','LineWidth',2,'MarkerFaceColor','w')
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

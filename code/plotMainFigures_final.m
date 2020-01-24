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


for i=1:size(compiledData_MVA,1)
    temp=nanmedian(compiledData_MVA(i).vds,3);
    temp2=compiledData_MVA(i).fooof_peak_Fr;
    between_factor(i+10)=1;
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            MVA_VDS(i,j,k)=temp(j,k);
            MVA_VEP(i,j,k)=temp2(j,k);
        end
    end
end

for i=1:size(compiledData_HAF,1)
    temp=nanmedian(compiledData_HAF(i).vds,3);
    temp2=compiledData_HAF(i).fooof_peak_Fr;
    between_factor(i+10)=1;
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            HAF_VDS(i,j,k)=temp(j,k);
            HAF_VEP(i,j,k)=temp2(j,k);     
        end
    end
end

%% Fit temporal frequency response function curves for bootstrap analysis


% % MVA VDS magno
% ydata=squeeze(MVA_VDS(:,1,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[0.5 4 1];
% [TTF_boot_mvaVDS_magno]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % MVA VDS parvo
% ydata=squeeze(MVA_VDS(:,2,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[2 2 1];
% [TTF_boot_mvaVDS_parvo]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % MVA VDS konio
% ydata=squeeze(MVA_VDS(:,3,[1:3 5]));
% xdata=[1.625 3.25 7.5 30];
% x0=[6 1 1];
% [TTF_boot_mvaVDS_konio]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % HAF VDS magno
% ydata=squeeze(HAF_VDS(:,1,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[0.5 4 1];
% [TTF_boot_hafVDS_magno]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % HAF VDS parvo
% ydata=squeeze(HAF_VDS(:,2,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[4 1 1];
% [TTF_boot_hafVDS_parvo]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % HAF VDS konio
% ydata=squeeze(HAF_VDS(:,3,[1:3 5]));
% xdata=[1.625 3.25 7.5 30];
% x0=[6 1 1];
% [TTF_boot_hafVDS_konio]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);


% 
% MVA VEP magno
ydata=squeeze(MVA_VEP(:,1,:));
xdata=[1.625 3.25 7.5 15 30];
x0=[0.5 4 1];
[TTF_boot_mvaVEP_magno]=bootstrp_ttf_fit(ydata,xdata,x0);

% % MVA VEP parvo
% ydata=squeeze(MVA_VEP(:,2,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[4 1 1];
% [TTF_boot_mvaVEP_parvo]=bootstrp_ttf_fit(ydata,xdata,x0);

% % MVA VEP konio
% ydata=squeeze(MVA_VEP(:,3,[1:3 5]));
% xdata=[1.625 3.25 7.5 30];
% x0=[5 2 1];
% [TTF_boot_mvaVEP_konio]=bootstrp_ttf_fit(ydata,xdata,x0);
% 
% HAF VEP magno
ydata=squeeze(HAF_VEP(:,1,:));
xdata=[1.625 3.25 7.5 15 30];
x0=[0.5 4 1];
[TTF_boot_hafVEP_magno]=bootstrp_ttf_fit(ydata,xdata,x0);

% % HAF VEP parvo
% ydata=squeeze(HAF_VEP(:,2,:));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[4 1 1];
% [TTF_boot_hafVEP_parvo]=bootstrp_ttf_fit(ydata,xdata,x0);

% % HAF VEP konio
% ydata=squeeze(HAF_VEP(:,3,[1:3 5]));
% xdata=[1.625 3.25 7.5 30];
% x0=[5 2 1];
% [TTF_boot_hafVEP_konio]=bootstrp_ttf_fit(ydata,xdata,x0);


% % VDS magno
% ydata=cat(1,squeeze(MVA_VDS(:,1,:)),squeeze(HAF_VDS(:,1,:)));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[0.5 4 1];
% [TTF_boot_VDS_magno]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % VDS parvo
% ydata=cat(1,squeeze(MVA_VDS(:,2,:)),squeeze(HAF_VDS(:,2,:)));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[2 2 1];
% [TTF_boot_VDS_parvo]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% % VDS konio
% ydata=cat(1,squeeze(MVA_VDS(:,3,[1:3 5])),squeeze(HAF_VDS(:,3,[1:3 5])));
% xdata=[1.625 3.25 7.5 30];
% x0=[6 1 1];
% [TTF_boot_VDS_konio]=bootstrp_ttf_fit_VDS(ydata,xdata,x0);
% 
% 
% VEP magno
ydata=cat(1,squeeze(MVA_VEP(:,1,:)),squeeze(HAF_VEP(:,1,:)));
xdata=[1.625 3.25 7.5 15 30];
x0=[0.5 4 1];
[TTF_boot_VEP_magno]=bootstrp_ttf_fit(ydata,xdata,x0);

% % VEP parvo
% ydata=cat(1,squeeze(MVA_VEP(:,2,:)),squeeze(HAF_VEP(:,2,:)));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[4 1 1];
% [TTF_boot_VEP_parvo]=bootstrp_ttf_fit(ydata,xdata,x0);

% % VEP konio
% ydata=cat(1,squeeze(MVA_VEP(:,3,[1:3 5])),squeeze(HAF_VEP(:,3,[1:3 5])));
% xdata=[1.625 3.25 7.5 30];
% x0=[5 2 1];
% [TTF_boot_VEP_konio]=bootstrp_ttf_fit(ydata,xdata,x0);


% figure(1)
% subplot(2,1,1)
% hold on
% errorbar(TTF_boot_mvaVDS_magno.peak_frequency(:,500),1.25,TTF_boot_mvaVDS_magno.peak_frequency(:,500)-TTF_boot_mvaVDS_magno.peak_frequency(:,50),TTF_boot_mvaVDS_magno.peak_frequency(:,950)-TTF_boot_mvaVDS_magno.peak_frequency(:,500),'ok','horizontal');
% errorbar(TTF_boot_mvaVDS_parvo.peak_frequency(:,500),1.25,TTF_boot_mvaVDS_parvo.peak_frequency(:,500)-TTF_boot_mvaVDS_parvo.peak_frequency(:,50),TTF_boot_mvaVDS_parvo.peak_frequency(:,950)-TTF_boot_mvaVDS_parvo.peak_frequency(:,500),'or','horizontal');
% errorbar(TTF_boot_mvaVDS_konio.peak_frequency(:,500),1.25,TTF_boot_mvaVDS_konio.peak_frequency(:,500)-TTF_boot_mvaVDS_konio.peak_frequency(:,50),TTF_boot_mvaVDS_konio.peak_frequency(:,950)-TTF_boot_mvaVDS_konio.peak_frequency(:,500),'ob','horizontal');
% errorbar(TTF_boot_hafVDS_magno.peak_frequency(:,500),1.75,TTF_boot_hafVDS_magno.peak_frequency(:,500)-TTF_boot_hafVDS_magno.peak_frequency(:,50),TTF_boot_hafVDS_magno.peak_frequency(:,950)-TTF_boot_hafVDS_magno.peak_frequency(:,500),'ok','horizontal');
% errorbar(TTF_boot_hafVDS_parvo.peak_frequency(:,500),1.75,TTF_boot_hafVDS_parvo.peak_frequency(:,500)-TTF_boot_hafVDS_parvo.peak_frequency(:,50),TTF_boot_hafVDS_parvo.peak_frequency(:,950)-TTF_boot_hafVDS_parvo.peak_frequency(:,500),'or','horizontal');
% errorbar(TTF_boot_hafVDS_konio.peak_frequency(:,500),1.75,TTF_boot_hafVDS_konio.peak_frequency(:,500)-TTF_boot_hafVDS_konio.peak_frequency(:,50),TTF_boot_hafVDS_konio.peak_frequency(:,950)-TTF_boot_hafVDS_konio.peak_frequency(:,500),'ob','horizontal');
% title(['Peak frequency - VDS'])
% ax=gca;ax.XScale='log';ax.XLim=[2 50];ax.YLim=[1 2];ax.TickDir='out';ax.Box='off';
% 
% subplot(2,1,2)
% hold on
% errorbar(TTF_boot_mvaVEP_magno.peak_frequency(:,500),1.25,TTF_boot_mvaVEP_magno.peak_frequency(:,500)-TTF_boot_mvaVEP_magno.peak_frequency(:,50),TTF_boot_mvaVEP_magno.peak_frequency(:,950)-TTF_boot_mvaVEP_magno.peak_frequency(:,500),'ok','horizontal');
% errorbar(TTF_boot_mvaVEP_parvo.peak_frequency(:,500),1.25,TTF_boot_mvaVEP_parvo.peak_frequency(:,500)-TTF_boot_mvaVEP_parvo.peak_frequency(:,50),TTF_boot_mvaVEP_parvo.peak_frequency(:,950)-TTF_boot_mvaVEP_parvo.peak_frequency(:,500),'or','horizontal');
% errorbar(TTF_boot_mvaVEP_konio.peak_frequency(:,500),1.25,TTF_boot_mvaVEP_konio.peak_frequency(:,500)-TTF_boot_mvaVEP_konio.peak_frequency(:,50),TTF_boot_mvaVEP_konio.peak_frequency(:,950)-TTF_boot_mvaVEP_konio.peak_frequency(:,500),'ob','horizontal');
% errorbar(TTF_boot_hafVEP_magno.peak_frequency(:,500),1.75,TTF_boot_hafVEP_magno.peak_frequency(:,500)-TTF_boot_hafVEP_magno.peak_frequency(:,50),TTF_boot_hafVEP_magno.peak_frequency(:,950)-TTF_boot_hafVEP_magno.peak_frequency(:,500),'ok','horizontal');
% errorbar(TTF_boot_hafVEP_parvo.peak_frequency(:,500),1.75,TTF_boot_hafVEP_parvo.peak_frequency(:,500)-TTF_boot_hafVEP_parvo.peak_frequency(:,50),TTF_boot_hafVEP_parvo.peak_frequency(:,950)-TTF_boot_hafVEP_parvo.peak_frequency(:,500),'or','horizontal');
% errorbar(TTF_boot_hafVEP_konio.peak_frequency(:,500),1.75,TTF_boot_hafVEP_konio.peak_frequency(:,500)-TTF_boot_hafVEP_konio.peak_frequency(:,50),TTF_boot_hafVEP_konio.peak_frequency(:,950)-TTF_boot_hafVEP_konio.peak_frequency(:,500),'ob','horizontal');
% title(['Peak frequency - VEP'])
% ax=gca;ax.XScale='log';ax.XLim=[2 50];ax.YLim=[1 2];ax.TickDir='out';ax.Box='off';

%% plot temporal frequency response functions (Figure 2)
VDS=calcVDS(compiledData_ALL,lb,ub);
VDS_mva=calcVDS(compiledData_MVA,lb,ub);
VDS_haf=calcVDS(compiledData_HAF,lb,ub);

[ttf_fitLMS_vds,TemporalFrequency_fitLMS_vds,paramsLMS_vds]=getTTFfits_VDS(VDS.LMS_M,TemporalFrequency,[0.5 4 1]);
[ttf_fitLM_vds,TemporalFrequency_fitLM_vds,paramsLM_vds]=getTTFfits_VDS(VDS.LM_M,TemporalFrequency,[2 2 1]);
[ttf_fitS_vds,TemporalFrequency_fitS_vds,paramsS_vds]=getTTFfits_VDS(VDS.S_M(1,[1:3 5]),TemporalFrequency([1:3 5]),[6 1 1]);

figure(2)
subplot(3,2,1)
hold on
markerline='ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LMS_M,VDS.LMS_CI,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS_vds,ttf_fitLMS_vds,'-k')
title(['Flicker discomfort MwA'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,3)
hold on
markerline='or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,VDS.LM_M,VDS.LM_CI,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM_vds,ttf_fitLM_vds,'-r')
ylabel(['FLicker dicomfort'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(3,2,5)
hold on
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),VDS.S_M(1,[1:3 5]),VDS.S_CI(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS_vds,ttf_fitS_vds,'-b')
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];


% Plot median psd for stimulus frequency
[LMSm, LMm, Sm, LMSci, LMci, Sci]=medianFooofFrequency(compiledData_ALL,lb,ub);

[ttf_fitLMS,TemporalFrequency_fitLMS,paramsLMS_vep]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM,paramsLM_vep]=getTTFfits(LMm,TemporalFrequency,[4 1 1]);
[ttf_fitS,TemporalFrequency_fitS,paramsS_vep]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

subplot(3,2,2)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='ok';markeredge=[0 0 0];markerface=[0 0 0];
plotWithErrorbars(TemporalFrequency,LMSm,LMSci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLMS,ttf_fitLMS,'-k')
title(['LMS'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,2,4)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='or';markeredge=[1 0 0];markerface=[1 0 0];
plotWithErrorbars(TemporalFrequency,LMm,LMci,markerline,markeredge,markerface)
plot(TemporalFrequency_fitLM,ttf_fitLM,'-r')
xlabel('Stimulus frequency')
title(['LM'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


subplot(3,2,6)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency([1:3 5]),Sm(1,[1:3 5]),Sci(:,[1:3 5]),markerline,markeredge,markerface)
plot(TemporalFrequency_fitS,ttf_fitS,'-b')
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];

% All Konio data including 15 Hz
figure(20)
subplot(2,1,1)
hold on
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency,VDS.S_M(1,:),VDS.S_CI(:,:),markerline,markeredge,markerface)
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[0 10];

subplot(2,1,2)
hold on
plot([0.95 35],[0 0],'--','Color',[0.8 0.8 0.8])
markerline='ob';markeredge=[0 0 1];markerface=[0 0 1];
plotWithErrorbars(TemporalFrequency,Sm(1,:),Sci(:,:),markerline,markeredge,markerface)
title(['S'])
ax=gca;ax.XScale='log';ax.XLim=[0.95 35];ax.YLim=[-0.002 0.022];


%% plot visual discomfort data as a function of VEP power at the stimulus
% frequency (FIgure 3)
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
ax=gca; ax.XLim=[0 0.011]; ax.YLim=[0 10]; ax.TickDir='out'; ax.Box='off';


clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci




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





function [TTF_boot]=bootstrp_ttf_fit(ydata,xdata,x0)
    
    Bootstat=bootstrp(1000,@nanmedian,ydata,1);
    for i=1:size(Bootstat,1)
        [ttf_fit,TF_fit,params]=getTTFfits(Bootstat(i,:),xdata,x0);
        temp=find(ttf_fit==nanmax(ttf_fit));
        TTF_boot.peak_frequency(i)=TF_fit(temp(1));
        TTF_boot.median_amplitude(i)=nanmedian(ttf_fit);
        TTF_boot.peak_amplitude(i)=nanmax(ttf_fit);
        TTF_boot.ttf_fit(i,:)=ttf_fit;
        TTF_boot.TF=TF_fit;
        TTF_boot.params(i,:)=params;
        
        figure(10)
        hold on
        plot(xdata,Bootstat(i,:),'og')
        plot(TF_fit,ttf_fit,'-g')
    end
    
        TTF_boot.peak_frequency=sort(TTF_boot.peak_frequency,2);
        TTF_boot.median_amplitude=sort(TTF_boot.median_amplitude,2);
        TTF_boot.peak_amplitude=sort(TTF_boot.peak_amplitude,2);
        TTF_boot.Bootstat=Bootstat;
        
%         pause
        clf
        
end


function [TTF_boot]=bootstrp_ttf_fit_VDS(ydata,xdata,x0)
    
    Bootstat=bootstrp(1000,@nanmedian,ydata,1);
    for i=1:size(Bootstat,1)
        [ttf_fit,TF_fit,params]=getTTFfits_VDS(Bootstat(i,:),xdata,x0);
        temp=find(ttf_fit==nanmax(ttf_fit));
        TTF_boot.peak_frequency(i)=TF_fit(temp(1));
        TTF_boot.median_amplitude(i)=nanmedian(ttf_fit);
        TTF_boot.peak_amplitude(i)=nanmax(ttf_fit);
        TTF_boot.ttf_fit(i,:)=ttf_fit;
        TTF_boot.TF=TF_fit;
        TTF_boot.params(i,:)=params;
        
        figure(10)
        hold on
        plot(xdata,Bootstat(i,:),'og')
        plot(TF_fit,ttf_fit,'-g')
    end
    
        TTF_boot.peak_frequency=sort(TTF_boot.peak_frequency,2);
        TTF_boot.median_amplitude=sort(TTF_boot.median_amplitude,2);
        TTF_boot.peak_amplitude=sort(TTF_boot.peak_amplitude,2);
        TTF_boot.Bootstat=Bootstat;
                
%         pause
        clf

        
end
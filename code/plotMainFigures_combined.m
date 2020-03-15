% Load compiled data across subjects
savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

filenameComp=fullfile(savePath,'allSubjects.mat');

open(filenameComp)
compiledData_MVA=ans.compiledData_MVA;
compiledData_HAF=ans.compiledData_HAF;
compiledData_ALL=ans.compiledData_HAF;
scoreTable_MVA=ans.scoreTable_MVA;
scoreTable_HAF=ans.scoreTable_HAF;
scoreTable=ans.scoreTable;

TemporalFrequency=[1.625 3.25 7.5 15 30];
lb=50;
ub=950;

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
    between_factor(i)=2;
    for j=1:size(temp,1)
        for k=1:size(temp,2)
               if j==3 && k==4
                wf_MVA_VDS(i,j,k)=temp(j,k);
                wf_MVA_VEP(i,j,k)=temp2(j,k);
                disp('ignore 15Hz S')
               else
                VDS_all=cat(2,VDS_all,temp(j,k));
                VEP_all=cat(2,VEP_all,temp2(j,k));
                group=cat(2,group,{'MVA'});
                TF=cat(2,TF,TemporalFrequency(k));
                
                wf_MVA_VDS(i,j,k)=temp(j,k);
                wf_MVA_VEP(i,j,k)=temp2(j,k);

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
end

for i=1:size(compiledData_HAF,1)
    temp=nanmedian(compiledData_HAF(i).vds,3);
    temp2=compiledData_HAF(i).fooof_peak_Fr;
    between_factor(i+10)=1;
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            if j==3 && k==4
                disp('ignore 15Hz S')
                wf_HAF_VDS(i,j,k)=temp(j,k);
                wf_HAF_VEP(i,j,k)=temp2(j,k);
                else
                    VDS_all=cat(2,VDS_all,temp(j,k));
                    VEP_all=cat(2,VEP_all,temp2(j,k));
                    group=cat(2,group,{'HAF'});
                    TF=cat(2,TF,TemporalFrequency(k));

                    wf_HAF_VDS(i,j,k)=temp(j,k);
                    wf_HAF_VEP(i,j,k)=temp2(j,k);
                    
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
end

figure(100)
temp=1:120;
temp2=121:length(VDS_all);
errorbar(1,median(VDS_all(temp)),std(VDS_all(temp))./sqrt(length(temp)-1),'ok')
hold on
errorbar(1,median(VDS_all(temp2)),std(VDS_all(temp2))./sqrt(length(temp2)-1),'ob')


%% Fit temporal frequency response function curves

% for averaged data for VDS
[ttf_fitLMS_vds,TemporalFrequency_fitLMS_vds,paramsLMS_vds]=getTTFfits_VDS(VDS.LMS_M,TemporalFrequency,[0.5 4 1]);
[ttf_fitLM_vds,TemporalFrequency_fitLM_vds,paramsLM_vds]=getTTFfits_VDS(VDS.LM_M,TemporalFrequency,[2 2 1]);
[ttf_fitS_vds,TemporalFrequency_fitS_vds,paramsS_vds]=getTTFfits_VDS(VDS.S_M(1,[1:3 5]),TemporalFrequency([1:3 5]),[6 1 1]);

% fit per subject
konio_TF=[1.625 3.25 7.5 30];
figure(300)
for x=1:size(compiledData_ALL,1)
    vds_subject_konio=squeeze(nanmedian(compiledData_ALL(x).vds(3,[1:3 5],:),3));
    vep_subject_konio=compiledData_ALL(x).fooof_peak_Fr(3,[1:3 5]);
    if sum(vds_subject_konio)==0
        vds_sim_15Hz(x,1)=0;
    else
        [ttf_fit_vds_sim,TemporalFrequency_fit_vds_sim,params_vds_sim]=getTTFfits_VDS(vds_subject_konio,konio_TF,[7 1 1]);
        sim_15Hz=find(TemporalFrequency_fit_vds_sim>15);
        vds_sim_15Hz(x,1)=ttf_fit_vds_sim(sim_15Hz(1));
        
        subplot(2,1,1)
        plot(konio_TF,vds_subject_konio,'ob')
        hold on
        plot(TemporalFrequency_fit_vds_sim,ttf_fit_vds_sim,'b')
        plot(15,ttf_fit_vds_sim(sim_15Hz(1)),'or')
        hold off
    end
    [ttf_fit_vep_sim,TemporalFrequency_fit_vep_sim,params_vep_sim]=getTTFfits(vep_subject_konio,konio_TF,[7 1 1]);
    sim_15Hz=find(TemporalFrequency_fit_vep_sim>15);
    vep_sim_15Hz(x,1)=ttf_fit_vep_sim(sim_15Hz(1));
    
    subplot(2,1,2)
    plot(konio_TF,vep_subject_konio,'ob')
    hold on
    plot(TemporalFrequency_fit_vep_sim,ttf_fit_vep_sim,'b')
    plot(15,ttf_fit_vep_sim(sim_15Hz(1)),'or')
%     pause
    hold off
   
end


% Boostrap analysis fit

% MwVA VDS
% ydata=squeeze(nanmedian(compiledData_MVA.vds(3,[1:3 5],:),3));
% xdata=[1.625 3.25 7.5 15 30];
% x0=[0.5 4 1];
% [TTF_boot]=bootstrp_ttf_fit(ydata,xdata,x0);


% %% Mixed effects repeated measures ANOVA
% withinfactor_VDS=cat(1,wf_MVA_VDS,wf_HAF_VDS);
% withinfactor_VEP=cat(1,wf_MVA_VEP,wf_HAF_VEP);
% between_factor=between_factor';
% 
% bf_names={'group'};
% wf_names={'direction','TF'};
% wf_variables=fullfact([3 5]);
% wf_tbl=array2table(wf_variables,'VariableNames',wf_names);
% wf_tbl.TF=categorical(wf_tbl.TF);
% wf_tbl.direction=categorical(wf_tbl.direction);
% 
% temp_vds=withinfactor_VDS(:,:);
% temp_vep=withinfactor_VEP(:,:);
% temp_vds_sim15=cat(2,temp_vds(:,1:11),vds_sim_15Hz,temp_vds(:,13:end));
% temp_vep_sim15=cat(2,temp_vep(:,1:11),vep_sim_15Hz,temp_vep(:,13:end));
% var_names=cellstr([{'y1'};{'y2'};{'y3'};{'y4'};{'y5'};{'y6'};{'y7'};{'y8'};{'y9'};{'y10'};{'y11'};{'y12'};{'y13'};{'y14'};{'y15'};{'group'}]);
% VDS_tbl=array2table([temp_vds between_factor],'VariableNames',var_names);
% VEP_tbl=array2table([temp_vep between_factor],'VariableNames',var_names);
% VDS_tbl_sim15=array2table([temp_vds_sim15 between_factor],'VariableNames',var_names);
% VEP_tbl_sim15=array2table([temp_vep_sim15 between_factor],'VariableNames',var_names);
% 
% % with real 15 Hz data
% rm_VDS=fitrm(VDS_tbl,'y1-y15~group','WithinDesign',wf_tbl);
% ranova_tbl_VDS_between=ranova(rm_VDS);
% ranova_tbl_VDS_within=ranova(rm_VDS,'WithinModel','direction+TF');
% 
% rm_VEP=fitrm(VEP_tbl,'y1-y15~group','WithinDesign',wf_tbl);
% ranova_tbl_VEP_between=ranova(rm_VEP);
% ranova_tbl_VEP_within=ranova(rm_VEP,'WithinModel','direction+TF');
% 
% % with simulated 15 Hz data
% rm_VDS_sim=fitrm(VDS_tbl_sim15,'y1-y15~group','WithinDesign',wf_tbl);
% ranova_tbl_VDS_between_sim15=ranova(rm_VDS_sim);
% ranova_tbl_VDS_within_sim15=ranova(rm_VDS_sim,'WithinModel','direction+TF');
% 
% rm_VEP_sim=fitrm(VEP_tbl_sim15,'y1-y15~group','WithinDesign',wf_tbl);
% ranova_tbl_VEP_between_sim15=ranova(rm_VEP_sim);
% ranova_tbl_VEP_within_sim15=ranova(rm_VEP_sim,'WithinModel','direction+TF');

%% plot temporal frequency response functions
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

[ttf_fitLMS,TemporalFrequency_fitLMS,paramsLMS_vep]=getTTFfits(LMSm,TemporalFrequency,[1 2 1]);
[ttf_fitLM,TemporalFrequency_fitLM,paramsLM_vep]=getTTFfits(LMm,TemporalFrequency,[4 1 1]);
[ttf_fitS,TemporalFrequency_fitS,paramsS_vep]=getTTFfits(Sm([1:3 5]),TemporalFrequency([1:3 5]),[6 2 1]);

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
ax=gca; ax.XLim=[0 0.011]; ax.YLim=[0 10]; ax.TickDir='out'; ax.Box='off';


clear LMSm LMm Sm BKGDm LMSci LMci Sci BKGDci


% plot luminance 30 Hz VEP response as a function of headache frequency
figure(6)
for i=1:size(compiledData_MVA,1)
    X_mva(i)=table2array(compiledData_MVA(i).subject(:,6));
    Y_mva(i)=compiledData_MVA(i).fooof_peak_Fr(1,5);
end

[Rmva_vepHAf,Pmva_vepHAf]=corrcoef(X_mva,Y_mva);
[b,stats]=robustfit(X_mva,Y_mva);

[LMSm, LMm, Sm, LMSci, LMci, Sci]=medianFooofFrequency(compiledData_HAF,lb,ub);
plot(X_mva,Y_mva,'ok','MarkerFaceColor','k','MarkerSize',8)
hold on
lsline
X=0:1:30;
plot(X,b(1)+b(2)*X,'--k')
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





function [TTF_boot]=bootstrp_ttf_fit(ydata,xdata,x0)
    
    Bootstat=bootstrp(1000,@nanmedian,ydata,1);
    for i=1:size(Bootstat,1)
        [ttf_fit,TF_fit,fit_params]=getTTFfits(bootstat,xdata,x0);
        temp=find(ttf_fit==nanmax(ttf_fit));
        TTF_boot.peak_frequency(i)=TF_fit(temp(1));
        TTF_boot.median_amplitude(i)=nanmedian(ttf_fit);
        TTF_boot.peak_amplitude(i)=nanmax(ttf_fit);
        TTF_boot.ttf_fit(i,:)=ttf_fit;
        TTF_boot.TF=TF_fit;
    end
    
        TTF_boot.peak_frequency=sort(TTF_boot.peak_frequency,1);
        TTF_boot.median_amplitude=sort(TTF_boot.median_amplitude,1);
        TTF_boot.peak_amplitude=sort(TTF_boot.peak_amplitude,1);
end
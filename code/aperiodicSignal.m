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
            MVA_nulling(counter_MVA,:)=compiledData.nulling;
        else if compiledData(1).group=='HA f'
                counter_HAF=counter_HAF+1;
               compiledData_HAF(counter_HAF,:)=compiledData;
               HAF_nulling(counter_HAF,:)=compiledData.nulling;
            else
                disp('error: group');
            end
        end

end

clear ans compiledData counter_HAF counter_MVA observerID x filenameComp

% Calculate anova of aperiodic signal
ap_all=[];
group=[];
PRP=[];
TF=[];
freq=[];
for i=1:size(compiledData_MVA,1)
    for j=1:3
        for k=1:5
            if j==3 && k==4
                disp('ignore 15Hz S')
            else
                temp=compiledData_MVA(i).fooof_results(j,k).bg_fit;
                ap_all=cat(2,ap_all,sum(temp));
                MVA_ap(i,j,k,:)=10.^temp;
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
end

for i=1:size(compiledData_HAF,1)
    for j=1:3
        for k=1:5
            if j==3 && k==4
                disp('ignore 15Hz S')
            else
                temp=compiledData_HAF(i).fooof_results(j,k).bg_fit;
                ap_all=cat(2,ap_all,sum(temp));
                HAF_ap(i,j,k,:)=10.^temp;
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
end

ap_anova=anovan(ap_all,{group,TF,PRP},'model','interaction','varnames',{'group','Temporal frequency','channel'});
for xx=1:size(MVA_ap,1)
    MVA_ap_subject(xx,:)=nanmedian(reshape(squeeze(MVA_ap(xx,:,:,:)),15,size(MVA_ap,4)));
end
for xx=1:size(HAF_ap,1)
    HAF_ap_subject(xx,:)=nanmedian(reshape(squeeze(HAF_ap(xx,:,:,:)),15,size(HAF_ap,4)));
end
freqs=compiledData_MVA(1).fooof_results(1,1).freqs;

figure(15)
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
ax.XLim=[0 60];
ax.YLim=[0 0.005];


% % plot median by subject
% subplot(1,3,1)
% hold on
% plot(freqs,MVA_ap_subject(xx,:),'k')
% title('Migraine with visual aura')
% ylabel('amplitude at stimulus frequency (mV)')
% ax=gca;
% ax.TickDir='out';
% ax.Box='off';
% ax.XLim=[0 35];
% ax.YLim=[-4 -2];
    
% % Plot nulling
% figure(10)
% hold on
% plot([1 2],MVA_nulling,'o','MarkerEdgeColor',[0.5 0.5 0.5])
% plot([1 2],median(MVA_nulling,1),'ok','MarkerFaceColor','k')
% plot([1.5 2.5],HAF_nulling,'o','MarkerEdgeColor',[0.5 0.5 0.5])
% plot([1.5 2.5],median(HAF_nulling,1),'ok','MarkerFaceColor','k')
% plot([0 1 3],[0 0 0],'--')
% [hLMn, pLMn]=ttest2(HAF_nulling(:,1),MVA_nulling(:,1));
% [hSn, pSn]=ttest2(HAF_nulling(:,2),MVA_nulling(:,2));
% title(['Nulling values, pLM=' num2str(pLMn) ', pS=' num2str(pSn)])
% ylabel('RGB adjustment values')
% ax=gca;
% ax.TickDir='out';
% ax.Box='off';
% ax.XTick=[1 1.5 2 2.5];
% ax.XTickLabel={'L-M MwA','L-M HAf','S MwA','S HAf'};
% ax.XLim=[0.5 3];
% ax.YLim=[-0.12 0.02];






% % Calculate anova of aperiodic signal
% ap_all=[];
% group=[];
% PRP=[];
% TF=[];
% freq=[];
% for i=1:size(compiledData_MVA,1)
%     for j=1:3
%         for k=1:5
%             temp=compiledData_MVA(i).fooof_results(j,k).bg_fit;
%             ap_all=cat(2,ap_all,temp);
%             for l=1:length(temp)
%                 group=cat(2,group,{'MVA'});
%                 TF=cat(2,TF,TemporalFrequency(k));
%                 freq=cat(2,freq,compiledData_MVA(i).fooof_results(j,k).freqs(l));
%         
%                 switch j
%                     case 1
%                         PRP=cat(2,PRP,{'LMS'});
%                     case 2
%                         PRP=cat(2,PRP,{'LM'});
%                     case 3
%                         PRP=cat(2,PRP,{'S'});
%                 end
%             end
%         end
%     end
% end
% 
% for i=1:size(compiledData_HAF,1)
%     for j=1:3
%         for k=1:5
%             temp=compiledData_HAF(i).fooof_results(j,k).bg_fit;
%             ap_all=cat(2,ap_all,temp);
%             for l=1:length(temp)
%                 group=cat(2,group,{'HAF'});
%                 TF=cat(2,TF,TemporalFrequency(k));
%                 freq=cat(2,freq,compiledData_HAF(i).fooof_results(j,k).freqs(l));
% 
%                 switch j
%                     case 1
%                         PRP=cat(2,PRP,{'LMS'});
%                     case 2
%                         PRP=cat(2,PRP,{'LM'});
%                     case 3
%                         PRP=cat(2,PRP,{'S'});
%                 end
%             end
%         end
%     end
% end
% 
% ap_anova=anovan(ap_all,{group,TF,PRP,freq},'model','interaction','varnames',{'group','Temporal frequency','channel','freq'});

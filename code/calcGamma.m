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

% % observer ID for each subject to be compiled (MwA only)
% subjects=[;'MELA_0131';'MELA_0175';...
%     'MELA_0170';'MELA_0174';...
%     'MELA_0194';'MELA_0179';...
%     'MELA_0209';'MELA_0207';...
%     'MELA_0213';'MELA_0208'];


TemporalFrequency=[1.625 3.25 7.5 15 30];
lb=50;
ub=950;
window=750;
Fs=2000;
C=['k';'r';'b'];
L=0.5:0.5:3;

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
        VEP_data(x,:,:,:)=squeeze(nanmean(compiledData.vep_Fr,3));

end

clear ans compiledData observerID x filenameComp temp*

figure(1)
for x=1:size(VEP_data,1)
    for y=1:size(VEP_data,2)
        for z=1:size(VEP_data,3)
            vep=squeeze(squeeze(VEP_data(x,y,z,:)));
            VEP(x,y,z,:)=vep;
%             [psd(x,y,z,:),freqs]=pwelch(vep',window,[],[],Fs);
            [psd(x,y,z,:),freqs]=pmtm(vep',4,length(vep),Fs);
        end
    end
end

clear x y z

vepM=squeeze(nanmean(VEP,1));
psdM=squeeze(nanmean(psd,1));

figure(2)
for x=1:size(psdM,1)
    for y=1:size(psdM,2)
        subplot(1,3,x)
        hold on
        plot(freqs,squeeze(psdM(x,y,:)),['-' C(x)],'LineWidth',L(y))
        ax=gca;ax.XLim=[0 100];ax.Box='off';ax.TickDir='out';ax.YScale='log';ax.YLim=[0.00000001 0.0001];
        if y==1
             plot([30 30],[0.000000001 0.001],'--','Color',[0.5 0.5 0.5])
             plot([80 80],[0.000000001 0.001],'--','Color',[0.5 0.5 0.5])
        end
    end
end

figure(3)
vep_plot=squeeze(vepM(1,3,:));
subplot(1,2,1)
plot((1:1:length(vep_plot))/Fs,vep_plot)
ax=gca;ax.Box='off';ax.TickDir='out';
subplot(1,2,2)
[S,F,T]=spectrogram(vep_plot',[],[],[],Fs);
helperCWTTimeFreqPlot(S,T,F,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.YLim=[0 100];
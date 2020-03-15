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

savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

load('/Users/carlynpattersongentile/Documents/MATLAB/Carlyn/vepMELA_subjectInfo.mat')

for x=1:size(subjects,1)
    
        observerID=subjects(x,:);
        
        filenameComp=fullfile(savePath,[observerID 'allChannels_1625w.mat']);
        open(filenameComp);
        compiledData_1625w=ans.compiledData;
        clear ans filenameComp
        
        filenameComp=fullfile(savePath,[observerID 'allChannels_325w.mat']);
        open(filenameComp);
        compiledData_325w=ans.compiledData;
        clear ans filenameComp
        
        filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);
        open(filenameComp);
        compiledData=ans.compiledData;
        clear ans filenameComp
        
        % replace 1.625Hz and 3.25Hz stimuli with appropriate windowing
        compiledData.fooof_peak_Fr(:,1)=compiledData_1625w.fooof_peak_Fr(:,1);
        compiledData.fooof_peak_harmonics(:,1)=compiledData_1625w.fooof_peak_harmonics(:,1);
        compiledData.fooof_peak_harmonics_freq(:,1)=compiledData_1625w.fooof_peak_harmonics_freq(:,1);
        compiledData.fooof_results(:,1)=compiledData_1625w.fooof_results(:,1);
        
        compiledData.fooof_peak_Fr(:,2)=compiledData_325w.fooof_peak_Fr(:,2);
        compiledData.fooof_peak_harmonics(:,2)=compiledData_325w.fooof_peak_harmonics(:,2);
        compiledData.fooof_peak_harmonics_freq(:,2)=compiledData_325w.fooof_peak_harmonics_freq(:,2);
        compiledData.fooof_results(:,2)=compiledData_325w.fooof_results(:,2);
        
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

clear compiledData counter_HAF counter_MVA observerID x filenameComp

filenameComp=fullfile(savePath,'allSubjects.mat');

save(filenameComp,'compiledData_MVA','scoreTable_MVA','compiledData_HAF','scoreTable_HAF','compiledData_ALL','scoreTable');
clear
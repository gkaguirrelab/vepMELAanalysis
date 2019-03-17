% Compile data across subjects

% observer ID for each subject to be compiled
subjects=['MELA_0121';'MELA_0131';...
    'MELA_0181';'MELA_0167';...
    'MELA_0187';'MELA_0175';...
    'MELA_0170';'MELA_0169'];

counter_MVA=0;
counter_HAF=0;
A=[1.625 3.25 7.5 15 30];

savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis','allChannels');

for x=1:size(subjects,1)
    
        observerID=subjects(x,:);
        savePath=fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
        'vepMELAanalysis','allChannels');

        filenameComp=fullfile(savePath,[observerID 'allChannels.mat']);

        open(filenameComp);
        compiledData=ans.compiledData;

        if compiledData(1).group=='MWVA'
            counter_MVA=counter_MVA+1;
            MVA_vds(counter_MVA,:,:,:)=compiledData.vds;
            MVA_fooof_fr(counter_MVA,:,:)=compiledData.fooof_peak_Fr;
        else if compiledData(1).group=='HA f'
                counter_HAF=counter_HAF+1;
                HAF_vds(counter_HAF,:,:,:)=compiledData.vds;
                HAF_fooof_fr(counter_HAF,:,:)=compiledData.fooof_peak_Fr;
            else
                disp('error: group');
            end
        end

end

    % Plot TFF by stimulus frequency
    
    
    figure(5)
    MVA_fooof_M=squeeze(mean(MVA_fooof_fr,1));
    MVA_fooof_std=squeeze(std(MVA_fooof_fr,[],1));
    HAF_fooof_M=squeeze(mean(HAF_fooof_fr,1));
    HAF_fooof_std=squeeze(std(HAF_fooof_fr,[],1));
    for y=1:size(MVA_vds,2)
        switch y
            case 1
                color='k';
                flicker_stim=1:5;
            case 2
                color='r';
                flicker_stim=1:5;
            case 3
                color='b';
                flicker_stim=[1:3 5];
        end
        subplot(2,2,1)
        errorbar(A(flicker_stim),MVA_fooof_M(y,flicker_stim),MVA_fooof_std(y,flicker_stim),['-o' color])
        hold on
        title(observerID)
        ylabel('power spectra for stimulus frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[0 0.15];
        
        subplot(2,2,2)
        errorbar(A(flicker_stim),HAF_fooof_M(y,flicker_stim),HAF_fooof_std(y,flicker_stim),['-.' color])
        hold on
        title(observerID)
        ylabel('power spectra for stimulus frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 35];
        ax.YLim=[0 0.15];
    end
 
    
    % Plot mean Visual discomfort data
    figure(5)
    MVA_vds=squeeze(nanmedian(MVA_vds,4));
    HAF_vds=squeeze(nanmedian(HAF_vds,4));
    for y=1:size(MVA_vds,2)
        switch y
            case 1
                color='k';
                flicker_stim=1:5;
            case 2
                color='r';
                flicker_stim=1:5;
            case 3
                color='b';
                flicker_stim=[1:3 5];
        end
        subplot(2,2,4)
        VDS_MVA=squeeze(nanmedian(MVA_vds,1));
        VDSstd_MVA=squeeze(nanstd(MVA_vds,[],1));
        errorbar(A(flicker_stim),VDS_MVA(y,flicker_stim),VDSstd_MVA(y,flicker_stim),['-o' color])
        hold on
        ylabel('visual discomfort scale')
        xlabel('temporal frequency of stimulus')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 65];
        ax.YLim=[0 11];
    
        subplot(2,2,3)
        VDS_HAF=squeeze(nanmedian(HAF_vds,1));
        VDSstd_HAF=squeeze(nanstd(HAF_vds,[],1));
        errorbar(A(flicker_stim),VDS_HAF(y,flicker_stim),VDSstd_HAF(y,flicker_stim),['-.' color])
        hold on
                ylabel('visual discomfort scale')
        xlabel('temporal frequency of stimulus')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XScale='log';
        ax.XLim=[0.95 65];
        ax.YLim=[0 11];
    end

    
% runFOOOF


% :size(PSD,1)
for x=1
    vep_FrM=squeeze(nanmean(processedVEPdata(x).vep_Fr,2));
    Fs=parsedVEPdata(x).Fs;
    for y=1
        % FOOOF inputs must be row vectors
        [psd, freqs] = pwelch(vep_FrM(y,:),1500,[],[],Fs);
        psd=psd'; freqs=freqs';

        % FOOOF settings
        f_range=[1, 150];
        settings=fooof_check_settings(settings);
        settings.background_mode='knee';

        % Run FOOOF
        fooof_results=fooof(freqs,psd,f_range,settings);

        model_fit=fooof_get_model(fooof_results);
    end
end
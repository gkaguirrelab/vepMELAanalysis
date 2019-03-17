function [fooof_results]=runFOOOF(vep_Fr,Fs,norm_vep)

switch length(size(vep_Fr))
    case 2
        vep_FrM=nanmean(vep_Fr,1)./norm_vep;
    case 3
        vep_FrM=squeeze(nanmean(vep_Fr,2))./norm_vep;
end

for y=1:size(vep_FrM,1)
    % FOOOF inputs must be row vectors
    [psd, freqs]=pwelch(vep_FrM(y,:),1500,[],[],Fs);
    
    if size(psd,1)>1
        psd=psd'; freqs=freqs';
    end

    % FOOOF settings
    f_range=[1, 40];
    settings=struct('peak_width_limits', [1.8, 12],'background_mode', 'fixed',...
        'min_peak_amplitude', 0.001,'peak_threshold', 1.0);

    % Run FOOOF
    fooof_results(y,:)=fooof(freqs,psd,f_range,settings,'return_model');
   
    
%     figure(20)
%     plot(fooof_results(y,:).freqs,fooof_results(y,:).power_spectrum,'-k')
%     hold on
%     plot(fooof_results(y,:).freqs,fooof_results(y,:).fooofed_spectrum,'-r')
%     title(num2str(y))
%     ax=gca;
%     ax.Box='off';
%     ax.TickDir='out';
%     ax.YLim=[-8 -0.1];
%     hold off
%     pause
end

end
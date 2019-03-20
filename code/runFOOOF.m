function [fooof_results]=runFOOOF(vep_Fr,Fs,window)


switch length(size(vep_Fr))
    case 2
        vep_FrM=nanmedian(vep_Fr,1);
        A=0;
    case 3
        vep_FrM=squeeze(nanmedian(vep_Fr,2));
        A=[1.625 3.25 7.5 15 30];
end

for y=1:size(vep_FrM,1)
    % FOOOF inputs must be row vectors
    [psd, freqs]=pwelch(vep_FrM(y,:),window,[],[],Fs);
    
    if size(psd,1)>1
        psd=psd'; freqs=freqs';
    end

    % FOOOF settings
    f_range=[1, 100];
    settings=struct('peak_width_limits', [2, 12],'background_mode', 'fixed',...
        'min_peak_amplitude', 0.5,'peak_threshold', 2,'max_n_peak', 10);

    % Run FOOOF
    fooof_results(y,:)=fooof(freqs,psd,f_range,settings,'return_model');
   
    
    figure(20)
    subplot(size(vep_FrM,1),1,y)
    plot(fooof_results(y,:).freqs,fooof_results(y,:).power_spectrum,'-k','LineWidth',2)
    hold on
    plot(fooof_results(y,:).freqs,fooof_results(y,:).fooofed_spectrum,'-r','LineWidth',2)
    plot(fooof_results(y,:).freqs,fooof_results(y,:).bg_fit,'--b','LineWidth',2)
    title(['stimulus frequency= ' num2str(A(y)) '. Error=' num2str(fooof_results(y,:).error) '. R squared=' num2str(fooof_results(y,:).r_squared)])
    ax=gca;
    ax.Box='off';
    ax.TickDir='out';
    ax.YLim=[-8 -1];
    hold off
end
pause

end
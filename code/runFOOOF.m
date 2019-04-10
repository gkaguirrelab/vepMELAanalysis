function [fooof_results,fooof_results5,fooof_results95]=runFOOOF(ttf_M,ttf_CI,Fs,dur_in_sec)

L=dur_in_sec*Fs;
freqs=Fs*(0:(L/2))/L;


for y=1:size(ttf_M,1)
    psd=ttf_M(y,:);
    
    switch size(ttf_M,1)
    case 1
        A=0;
        psd5=squeeze(ttf_CI(:,1))';
        psd95=squeeze(ttf_CI(:,2))';
    case 5
        A=[1.625 3.25 7.5 15 30];
        psd5=squeeze(ttf_CI(y,:,1));
        psd95=squeeze(ttf_CI(y,:,2));
    end

    % FOOOF settings
    f_range=[1, 100];
    settings=struct('peak_width_limits', [1.5, 12],'background_mode', 'knee',...
        'min_peak_amplitude', 0,'peak_threshold', 2);

    % Run FOOOF
    fooof_results(y,:)=fooof(freqs,psd,f_range,settings,'return_model');
    
    fooof_results5(y,:)=fooof(freqs,psd5,f_range,settings,'return_model');
    fooof_results95(y,:)=fooof(freqs,psd95,f_range,settings,'return_model');
   
    
%     figure(20)
%     subplot(size(ttf_M,1),1,y)
%     plot(fooof_results(y,:).freqs,fooof_results(y,:).power_spectrum,'-k','LineWidth',2)
%     hold on
%     plot(fooof_results(y,:).freqs,fooof_results(y,:).fooofed_spectrum,'-r','LineWidth',2)
%     plot(fooof_results(y,:).freqs,fooof_results(y,:).bg_fit,'--b','LineWidth',2)
%     title(['stimulus frequency= ' num2str(A(y)) '. Error=' num2str(fooof_results(y,:).error) '. R squared=' num2str(fooof_results(y,:).r_squared)])
%     ax=gca;
%     ax.Box='off';
%     ax.TickDir='out';
%     ax.YLim=[-5 -1];
%     hold off
end
% pause

end
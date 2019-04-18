function [fooof_results,fooof_results5,fooof_results95]=runFOOOF(ttf_M,ttf_CI,varargin)

%% Parse input
p = inputParser;
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30],@isnumeric);
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('dur_in_sec',1.5,@isnumeric);
p.addParameter('plot_all',false,@islogical);
p.addParameter('f_range',[1 100],@isnumeric);

p.parse(varargin{:});

L=p.Results.dur_in_sec*p.Results.Fs;
freqs=p.Results.Fs*(0:(L/2))/L;


for y=1:size(ttf_M,1)
    psd=ttf_M(y,:);
    
    switch size(ttf_M,1)
    case 1
        psd5=squeeze(ttf_CI(:,1))';
        psd95=squeeze(ttf_CI(:,2))';
    case 5
        psd5=squeeze(ttf_CI(y,:,1));
        psd95=squeeze(ttf_CI(y,:,2));
    end

    % FOOOF settings
    settings=struct('peak_width_limits', [1.5, 12],'background_mode', 'knee',...
        'min_peak_amplitude', 0,'peak_threshold', 2);

    % Run FOOOF
    fooof_results(y,:)=fooof(freqs,psd,p.Results.f_range,settings,'return_model');
    
    fooof_results5(y,:)=fooof(freqs,psd5,p.Results.f_range,settings,'return_model');
    fooof_results95(y,:)=fooof(freqs,psd95,p.Results.f_range,settings,'return_model');
   

    if p.Results.plot_all==1
        figure(20)
        subplot(size(ttf_M,1),1,y)
        plot(fooof_results(y,:).freqs,fooof_results(y,:).power_spectrum,'-k','LineWidth',2)
        hold on
        plot(fooof_results(y,:).freqs,fooof_results(y,:).fooofed_spectrum,'-r','LineWidth',2)
        plot(fooof_results(y,:).freqs,fooof_results(y,:).bg_fit,'--b','LineWidth',2)
        title(['stimulus frequency= ' num2str(p.Results.TemporalFrequency(y)) '. Error=' num2str(fooof_results(y,:).error) '. R squared=' num2str(fooof_results(y,:).r_squared)])
        ax=gca;
        ax.Box='off';
        ax.TickDir='out';
        ax.YLim=[-5 -1];
        hold off
    end
end

    if p.Results.plot_all==1
        pause
    end

end
function [fooof_results]=runFOOOF_bkgd(ttf_M,varargin)
%% Parse input
p = inputParser;
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('dur_in_sec',1.5,@isnumeric);
p.addParameter('f_range',[1 100],@isnumeric);

p.parse(varargin{:});

L=p.Results.dur_in_sec*p.Results.Fs;
freqs=p.Results.Fs*(0:(L/2))/L;


    for y=1:size(ttf_M,1)
        psd=ttf_M(y,:);
        
        % FOOOF settings
        settings=struct('peak_width_limits', [1.5, 12],'background_mode', 'knee',...
            'min_peak_amplitude', 0,'peak_threshold', 2);

        % Run FOOOF
        fooof_results(y,:)=fooof(freqs,psd,p.Results.f_range,settings,'return_model');

    end
end
function [ttf_VEP]=calcVEPttf(vep_Fr,varargin)
% A function that calculates Temporal transfer functions
% for Metropsis/VEP stimuli for individual observers.
%
% Syntax:
%  [ttf_VEP]=calcVEPttf(vep_Fr,Fs,dur_in_freq,varargin)
%
% Description:
%	This function accesses other functions to run all the analysis code for Metropsis/VEP stimuli.
%
%
% Output: saves structure VEP
%   ttf                   - Temporal transfer function

%% Parse input
p = inputParser;
p.addParameter('TemporalFreqency',[1.625 3.25 7.5 15 30 60],@isnumeric);
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('dur_in_freq',4000,@isnumeric);
p.addParameter('normalize',false,@islogical);

p.parse(varargin{:});

f=p.Results.Fs*(0:(p.Results.dur_in_freq/2))/p.Results.dur_in_freq; 

% Calculate fourier transform w is temporal frequency of the stimuli, x is
% the repeats (concatonated across sessions
for w=1:size(vep_Fr,1)
    for x=1:size(vep_Fr,2)
       % Fourier transform
        ft=fft(squeeze(vep_Fr(w,x,:)));
        P = abs(ft/p.Results.dur_in_freq);
        P_data(w,x,:) = P(1:p.Results.dur_in_freq/2+1);

        % select power spectra for each frequency
        temp=find(f>=p.Results.TemporalFreqency(x));
        P_dataFr(w,x)=squeeze(P_data(w,x,temp(1)));
    end

    subplot(1,2,1)
    plot(XX,squeeze(squeeze(mean(vep_Fr(w,x,:),3))))
    title(['frequency=' num2str(TF(x))]);
    xlabel('Time(s)')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.YLim=[-0.1 0.1];
    ax.XLim=[0 dur_in_sec];

    subplot(1,2,2)
    plot(f,squeeze(squeeze(mean(P_data(w,x,:,:),3))),'-k')
    hold on
    plot(TF(x),P_dataFr(w,x),'ob')
    ylabel('power spectra')
    xlabel('frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 130];
    ax.YLim=[0 0.02];
    pause
    hold off

end


if ('normalize')==1
    norm_ft=fft(mean(mean(vep_Fr,1),2));
    P=abs(norm_ft/p.Results.dur_in_freq);
    norm_ttf=P(1:p.Results.dur_in_freq/2+1);
    vep_Fr=vep_Fr-norm_ttf;
end

VEP_FrM=squeeze(mean(vep_Fr,2));    


for xx=1:size(VEP_Fr,2)
    ft=fft(VEP_FrM(xx,:));
    P=abs(ft/p.Results.dur_in_freq);
    P2=P(1:p.Results.dur_in_freq/2+1);

    figure(3)
    subplot(1,2,1)
    plot(XX,VEP_FrM(xx,:),'-k')
    title(['frequency=' num2str(TF(xx))]);
    xlabel('Time(s)')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.YLim=[-0.1 0.1];
    ax.XLim=[0 dur_in_sec];

    subplot(1,2,2)
    plot(f,P2,'-k')
    ylabel('power spectra')
    xlabel('frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 130];
    ax.YLim=[0 0.02];
    pause
    hold off
end


ttf_VEP.ttf=ttf;
ttf_VEP.ttfFr=ttfFr;
ttf_VEP.f=f;
end

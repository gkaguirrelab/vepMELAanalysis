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
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30],@isnumeric);
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('dur_in_sec',1.5,@isnumeric);
p.addParameter('plot_all',false,@islogical);

p.parse(varargin{:});

dur_in_freq=p.Results.dur_in_sec*p.Results.Fs;
f=p.Results.Fs*(0:(dur_in_freq/2))/dur_in_freq;
XX=(1:length(vep_Fr))/p.Results.Fs;
vep_FrM=squeeze(nanmean(vep_Fr,2));

for xx=1:size(vep_FrM,1)
    ft=fft(vep_FrM(xx,:));
    P=abs(ft/dur_in_freq);
    Bootstat=bootstrp(1000,@nanmean,squeeze(vep_Fr(xx,:,:)),1);
    for yy=1:size(Bootstat,1)
        boot_ft=fft(Bootstat(yy,:));
        P_boot=abs(boot_ft/dur_in_freq);
        ttf_boot=P_boot(1:dur_in_freq/2+1);
        temp=find(f>=p.Results.TemporalFrequency(xx));
        ttf_Fr_boot(:,yy)=max(ttf_boot(:,temp(1)-1:temp(1)));
    end
    
    ttf_Fr_boot=sort(ttf_Fr_boot);
    ttf_CI(xx,:)=ttf_Fr_boot(:,[50 950]);
    ttf_M(xx,:)=P(1:dur_in_freq/2+1);
    temp=find(f>=p.Results.TemporalFrequency(xx));
    ttf_FrM(xx,:)=max(ttf_M(xx,temp(1)-1:temp(1)));
    if p.Results.plot_all==1
        figure(3)
        subplot(1,2,1)
        plot(XX,vep_FrM(xx,:),'-k')
        title(['frequency=' num2str(p.Results.TemporalFrequency(xx))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 dur_in_freq/p.Results.Fs];

        subplot(1,2,2)
        plot(f,ttf_M(xx,:),'-k')
        hold on
        plot(p.Results.TemporalFrequency(xx),ttf_FrM(xx,:),'ob')
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
end


ttf_VEP.f=f;
ttf_VEP.ttf_M=ttf_M;
ttf_VEP.ttf_FrM=ttf_FrM;
ttf_VEP.ttf_CI=ttf_CI;
end

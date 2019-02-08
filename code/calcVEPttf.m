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
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30 60],@isnumeric);
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('dur_in_freq',4000,@isnumeric);
p.addParameter('normalize',false,@islogical);
p.addParameter('plot_all',false,@islogical);

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
        temp=find(f>=p.Results.TemporalFrequency(w));
        P_dataFr(w,x)=squeeze(P_data(w,x,temp(1,1)));
        clear temp
    end
    
    XX=(1:length(vep_Fr))/p.Results.Fs;
    
%     subplot(1,2,1)
%     plot(XX,squeeze(squeeze(mean(vep_Fr(w,x,:),3))))
%     title(['frequency=' num2str(p.Results.TemporalFrequency(w))]);
%     xlabel('Time(s)')
%     ax=gca;
%     ax.TickDir='out';
%     ax.Box='off';
%     ax.YLim=[-0.1 0.1];
%     ax.XLim=[0 p.Results.dur_in_freq/p.Results.Fs];
% 
%     subplot(1,2,2)
%     plot(f,squeeze(squeeze(mean(P_data(w,x,:,:),3))),'-k')
%     hold on
%     plot(p.Results.TemporalFrequency(w),P_dataFr(w,x),'ob')
%     ylabel('power spectra')
%     xlabel('frequency')
%     ax=gca;
%     ax.TickDir='out';
%     ax.Box='off';
%     ax.XLim=[0 130];
%     ax.YLim=[0 0.02];
%     pause
%     hold off

end

if p.Results.normalize==1
    clear P_dataFr
    norm_ft=fft(mean(mean(vep_Fr,1),2));
    norm_vep=mean(mean(vep_Fr,1),2);
    P=abs(norm_ft/p.Results.dur_in_freq);
    norm_ttf=P(1:p.Results.dur_in_freq/2+1);
    for x=1:size(P_data,1)
        for y=1:size(P_data,2)
            P_data(x,y,:)=P_data(x,y,:)-norm_ttf;
            temp=find(f>=p.Results.TemporalFrequency(x));
            P_dataFr(x,y)=squeeze(P_data(x,y,temp(1,1)));
            vep_Fr(x,y,:)=vep_Fr(x,y,:)-norm_vep;
        end
    end
    
    norm_ttf=squeeze(norm_ttf)';
    
    figure(10)
    plot(f,norm_ttf)
    ylabel('power spectra')
    xlabel('frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 130];
    ax.YLim=[0 0.02];
end

VEP_FrM=squeeze(mean(vep_Fr,2));    


for xx=1:size(vep_Fr,1)
    ft=fft(VEP_FrM(xx,:));
%     [Bootstat,Bootsam]=bootstrp(1000,@mean,parsed_VEP);
%     CI=sort(Bootstat);
%     CI95=Bootstat(950);
%     CI
    P=abs(ft/p.Results.dur_in_freq);
    ttf_M(xx,:)=P(1:p.Results.dur_in_freq/2+1);
    temp=find(f>=p.Results.TemporalFrequency(xx));
    ttf_FrM(xx,:)=ttf_M(xx,temp(1));
    if p.Results.plot_all==1
        figure(3)
        subplot(1,2,1)
        plot(XX,VEP_FrM(xx,:),'-k')
        title(['frequency=' num2str(p.Results.TemporalFrequency(xx))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 p.Results.dur_in_freq/p.Results.Fs];

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


ttf_VEP.ttf=P_data;
ttf_VEP.ttfFr=P_dataFr;
ttf_VEP.f=f;
ttf_VEP.vep_Fr=vep_Fr;
ttf_VEP.VEP_FrM=VEP_FrM;
ttf_VEP.ttf_M=ttf_M;
ttf_VEP.ttf_FrM=ttf_FrM;
end

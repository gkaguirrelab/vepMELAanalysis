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

L=p.Results.dur_in_sec*p.Results.Fs;
XX=(1:length(vep_Fr))/p.Results.Fs;
vep_FrM=squeeze(nanmedian(vep_Fr,2));
f=p.Results.Fs*(0:(L/2))/L;
f=f';
counter=1;

for xx=1:size(vep_FrM,1)
    psd_temp=fft(vep_FrM(xx,:));
    psd_temp=abs(psd_temp/L);
    ttf_M(xx,:)=psd_temp(:,1:(L/2)+1);
    Bootstat=bootstrp(100,@nanmedian,squeeze(vep_Fr(xx,:,:)),1);
    for yy=1:size(Bootstat,1)
        ttf_boot=fft(Bootstat(yy,:));
        ttf_boot=abs(ttf_boot/L);
        ttf_boot=ttf_boot(1:L/2+1);
        
        TTF_boot(xx,:,yy)=ttf_boot;
        temp=abs(f-p.Results.TemporalFrequency(xx));
        temp2=find(temp==min(temp));
        ttf_Fr_boot(:,yy)=max(ttf_boot(:,temp2));
    end
   
    TTF_boot=sort(TTF_boot,3);
    ttf_CI(xx,:,:)=TTF_boot(xx,:,[5 95]);
    
    temp=abs(f-p.Results.TemporalFrequency(xx));
    temp2=find(temp==min(temp));
    ttf_FrM(xx,:)=max(ttf_M(xx,temp2));
    ttf_Fr_boot=sort(ttf_Fr_boot);
    ttf_FrCI(xx,:)=ttf_Fr_boot(:,[5 95]);
    
   
    if p.Results.plot_all==1
        figure(3)
        subplot(5,2,counter)
        plot(XX,vep_FrM(xx,:),'-k')
        title(['frequency=' num2str(p.Results.TemporalFrequency(xx))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 p.Results.dur_in_sec];

        subplot(5,2,counter+1)
        plot(f,ttf_M(xx,:),'-k')
        hold on
        plot(p.Results.TemporalFrequency(xx),ttf_FrM(xx,:),'ob')
        ylabel('power spectra')
        xlabel('frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XLim=[0 100];
        ax.YLim=[0 0.01];
        pause
        hold off
        counter=counter+2;
    end
end

f=f';
ttf_VEP.f=f;
ttf_VEP.ttf_M=ttf_M;
ttf_VEP.ttf_CI=ttf_CI;
ttf_VEP.ttf_FrM=ttf_FrM;
ttf_VEP.ttf_FrCI=ttf_FrCI;
end

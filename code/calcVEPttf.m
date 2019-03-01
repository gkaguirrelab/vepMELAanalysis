function [ttf_VEP]=calcVEPttf(vep_Fr,vep_bkgd,varargin)
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
p.addParameter('dur_in_freq',4000,@isnumeric);
p.addParameter('normalize',false,@islogical);
p.addParameter('plot_all',false,@islogical);

p.parse(varargin{:});

f=p.Results.Fs*(0:(p.Results.dur_in_freq/2))/p.Results.dur_in_freq;
XX=(1:length(vep_Fr))/p.Results.Fs;

% Calculate fourier transform w is temporal frequency of the stimuli, x is
% the repeats (concatonated across sessions)
for w=1:size(vep_Fr,1)
    for x=1:size(vep_Fr,2)
       % Fourier transform
        ft=fft(squeeze(vep_Fr(w,x,:)));
        P = abs(ft/p.Results.dur_in_freq);
        P_data(w,x,:) = P(1:p.Results.dur_in_freq/2+1);
        
        % Plot each trial and exclude trials with noisy VEP signal
        figure(15)
        subplot(1,2,1)
        plot(XX,squeeze(vep_Fr(w,x,:)))
        title(['Frequency: ' num2str(p.Results.TemporalFrequency(w)) ' , Trial: ' num2str(x)])
        ax=gca;
        ax.YLim=[-0.5 0.5];
        ax.XLim=[0 p.Results.dur_in_freq/p.Results.Fs];
        ax.TickDir='out';
        ax.Box='off';
        hold off
        
        subplot(1,2,2)
        plot(f,squeeze(P_data(w,x,:)),'k')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        
        if max(diff(squeeze(vep_Fr(w,x,:))))>0.08 | max(diff(squeeze(vep_Fr(w,x,:))))<0.02
            title('bad')
        end

        hold off
        
        pause
        
        % select power spectra for each frequency
        temp=find(f>=p.Results.TemporalFrequency(w));
        P_dataFr(w,x)=max(squeeze(P_data(w,x,temp(1,1)-1:temp(1,1))))';
        clear temp
    end
  

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
end

VEP_FrM=squeeze(mean(vep_Fr,2));    


for xx=1:size(vep_Fr,1)
    ft=fft(VEP_FrM(xx,:));
    P=abs(ft/p.Results.dur_in_freq);
    Bootstat=bootstrp(1000,@mean,squeeze(vep_Fr(xx,:,:)),1);
    for yy=1:size(Bootstat,1)
        boot_ft=fft(Bootstat(yy,:));
        P_boot=abs(boot_ft/p.Results.dur_in_freq);
        ttf_boot=P_boot(1:p.Results.dur_in_freq/2+1);
        temp=find(f>=p.Results.TemporalFrequency(xx));
        ttf_Fr_boot(:,yy)=max(ttf_boot(:,temp(1)-1:temp(1)));
    end
    
    ttf_Fr_boot=sort(ttf_Fr_boot);
    ttf_CI(xx,:)=ttf_Fr_boot(:,[50 950]);
    ttf_M(xx,:)=P(1:p.Results.dur_in_freq/2+1);
    temp=find(f>=p.Results.TemporalFrequency(xx));
    ttf_FrM(xx,:)=max(ttf_M(xx,temp(1)-1:temp(1)));
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
ttf_VEP.ttf_CI=ttf_CI;
end

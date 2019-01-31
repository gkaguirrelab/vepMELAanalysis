function []=calcVEPttf(vep_Fr,vds_Fr,TF_trials)

A=unique(TF_trials);
f=Fs*(0:(dur_in_freq/2))/dur_in_freq; 

% Calculate fourier transform
for w=1:size(VDS_Fr,1)
    for x=1:size(VDS_Fr,2)
        for y=1:size(VDS_Fr,3)
           % Fourier transform
            ft=fft(squeeze(VEP_Fr(w,x,y,:)));
            P = abs(ft/dur_in_freq);
            P_data(w,x,y,:) = P(1:dur_in_freq/2+1);

            % select power spectra for each frequency
            temp=find(f>=A(x) & f<A(x)+diff(f(1:2)));
            P_dataFr(w,x,y)=squeeze(P_data(w,x,yy,temp(1)));
        end

        subplot(1,2,1)
        plot(XX,squeeze(squeeze(mean(VEP_Fr(w,x,:,:),3))))
        title(['frequency=' num2str(A(x))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 dur_in_sec];

        subplot(1,2,2)
        plot(f,squeeze(squeeze(mean(P_data(w,x,:,:),3))),'-k')
        hold on
        plot(60,squeeze(squeeze(mean(Noise60(w,x,:),3))),'or')
        ylabel('power spectra')
        xlabel('frequency')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XLim=[0 130];
        ax.YLim=[0 0.02];
        pause
        hold off

        yy=1;
    end
end
clear w yy x P ft

% collect 60Hz data from all freq to subtract out noise *** change this to
% normalize all data by averaging across frequencies
temp=find(f>=A(6) & f<A(6)+diff(f(1:2)));
Noise60(w,x,:)=squeeze(P_data(w,x,:,temp));

% Average data across repeats from all sessions
vep_Fr=[];
vds_Fr=[];
p_dataFr=[];
noise60=[];
for yy=1:size(VEP_Fr,1)
    vep_Fr=cat(2,vep_Fr,squeeze(VEP_Fr(yy,:,:,:)));
    vds_Fr=cat(2,vds_Fr,squeeze(VDS_Fr(yy,:,:)));
    p_dataFr=cat(2,p_dataFr,squeeze(P_dataFr(yy,:,:)));
    noise60=cat(2,noise60,Noise60(yy,:));
end
VEP_FrM=squeeze(mean(vep_Fr,2));


for xx=1:size(VEP_Fr,2)
    ft=fft(VEP_FrM(xx,:));
    P=abs(ft/dur_in_freq);
    P2=P(1:dur_in_freq/2+1);

    figure(3)
    subplot(1,2,1)
    plot(XX,VEP_FrM(xx,:),'-k')
    title(['frequency=' num2str(A(xx))]);
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

% collect 60Hz data from all freq to subtract out noise
temp=find(f>=A(6) & f<A(6)+diff(f(1:2)));
Noise60(w,x,:)=squeeze(P_data(w,x,:,temp));
        
figure(4)
% Plot TTF
subplot(2,1,1)
p_dataFr=squeeze(p_dataFr);
P_Fm=mean(p_dataFr,2);
P_Fm2=mean(p_dataFr(6,:),2)-mean(noise60);
P_Fstd=std(p_dataFr,[],2);
errorbar(A,cat(2,P_Fm(1:5)',P_Fm2),P_Fstd,'-ok')
hold on
errorbar(A(6),P_Fm(6),P_Fstd(6),'-ob')
title(expID)
ylabel('power spectra for stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 65];
ax.YLim=[0 0.02];

% Plot mean Visual discomfort data
subplot(2,1,2)
VDSm=nanmean(vds_Fr,2);
VDSstd=nanstd(vds_Fr,[],2);
errorbar(A,VDSm,VDSstd,'-ok')
ylabel('visual discomfort scale')
xlabel('temporal frequency of stimulus')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 65];
ax.YLim=[0 10];

% Plot TTF
figure(5)
hold on
p_dataFr=squeeze(p_dataFr);
P_Fm=mean(p_dataFr,2);
P_Fm2=mean(p_dataFr(6,:),2)-mean(noise60);
P_Fstd=std(p_dataFr,[],2);
errorbar(A,cat(2,P_Fm(1:5)',P_Fm2),P_Fstd,'-ob')
% hold on
% errorbar(A(6),P_Fm(6),P_Fstd(6),'-ob')
title(expID)
ylabel('power spectra for stimulus frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
ax.XLim=[0.95 65];
ax.YLim=[0 0.02];

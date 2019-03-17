% simulatedPSD

% Produce a 120Hz sampled stimulus to match the monitor
dur_in_sec=1.5; % duration of stimulus collection
Fs_stim=120; % sampling rate of the monitor
Fs=2000; % sampling rate of the brain
dur_in_freq=dur_in_sec*Fs_stim; % duration in sampling frequency of the monitor
sampling_monitor=0:1/Fs_stim:dur_in_sec;
sampling_VEP=0:1/Fs:dur_in_sec;
stim_freq=[1.625 3.25 7.5 15 30];
plot_color=['k','b','r','g','c','m','y'];
window=[750 1500];

for x=1:length(stim_freq)
    flicker=sin(2*pi*sampling_monitor*stim_freq(x));
    for y=1:length(window)
        for z=1:length(sampling_VEP)
            temp1=sampling_monitor-sampling_VEP(z);
            temp2=min(abs(temp1));
            temp3=find(abs(temp1)==temp2);
            simVEP(:,z)=flicker(:,temp3(1));
        end
        if y==1
            figure(1)
            plot(sampling_VEP,simVEP,'-k')
            hold on
            plot(sampling_monitor,flicker,'-r')
            title(['Frequency= ' num2str(stim_freq(x))])
            ax=gca;
            ax.Box='off';
            ax.TickDir='out';
            pause
            hold off
        end

        [psd, freqs] = pwelch(simVEP,window(y),[],[],Fs);
        figure(2)
        subplot(1,length(window),y)
        plot(freqs,psd,plot_color(x))
        hold on
        title(['Window= ' num2str(window(y))])
        ax=gca;
        ax.Box='off';
        ax.TickDir='out';
%         ax.XScale='log';
        ax.XLim=[0.1 152];
        ax.YLim=[0 1];
        pause
    end
end
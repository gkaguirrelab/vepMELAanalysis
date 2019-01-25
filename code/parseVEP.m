% parseVEP
% Code to organize VEP signal and parse trials based on input, and organize
% stimuli by stimulus frequency
clear; clc;

% load saved raw data
load 'test_VEP_data2.mat'

% VEP data is p.response(2,:)
% TTL pulse is p.response(1,:);
% timebase is p.timebase

% Bandpass filter for VEP signal
Fs=p.params.frequencyInHz;
lo=0.5; % low cut off frequency
hi=500; % high cut off frequency

d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',lo,...
    'HalfPowerFrequency2',hi,'SampleRate',Fs);

% Get rid of 60hz powerline hum
d2=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59,...
    'HalfPowerFrequency2',61,'SampleRate',Fs);
    
VEP_data=p.response(2,:);

VEP_data=filter(d,VEP_data);

VEP_data=filter(d2,VEP_data);

% VEP_data=p.response(2,:);

clear lo hi

% Find timestamp of TTL pulses
TTL=p.response(1,:);
D=100;
timestamp=p.timebase;
y=0;

for x=2:length(TTL)
    if TTL(1,x)>4 && TTL(1,(x-1))<4
        y=y+1;
        sync_pulse{y}=timestamp(1,x-D);
        sync_loc{y}=x-D;
        TTL_check{y}=TTL(1,x-D);
    end
end

sync_pulse=cell2mat(sync_pulse);
sync_loc=cell2mat(sync_loc);
TTL_check=cell2mat(TTL_check);

clear x y

% parse VEP data
dur_in_sec=2.5; % the length of the stimulus presentation
dur_in_freq=dur_in_sec*Fs;
repeat=7;

for x=1:length(sync_loc)
    parsed_VEP(x,:)=VEP_data(1,sync_loc(x):sync_loc(x)+dur_in_freq);
end

clear x

XX=(1:length(parsed_VEP))/Fs;
YY=mean(parsed_VEP,1);

figure(1)
plot(XX,YY,'-k')
title('Average across frequencies');
xlabel('Time(s)')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.YLim=[-0.1 0.1];
ax.XLim=[0 dur_in_sec];
hold off
% pause
% 
% for x=1:(size(parsed_VEP,1))
%     plot(XX,parsed_VEP(x,:))
%     hold off
%     pause
% end
% 
% clear x

% parse data by frequency
figure(2)

A=unique(trials);
VEP_Fr=zeros(length(A),repeat,length(parsed_VEP));
xx=1;
yy=1;
t=(0:dur_in_freq-1)*1/Fs;       % Time vector
f = Fs*(0:(dur_in_freq/2))/dur_in_freq; 


for x=1:length(A)
    for y=1:length(trials)
        if A(x)==trials(y)
           VEP_Fr(xx,yy,:)=parsed_VEP(y,:);
            ft=fft(parsed_VEP(y,:));
            P = abs(ft/dur_in_freq);
            P_data(xx,yy,:) = P(1:dur_in_freq/2+1);
           yy=yy+1;
        end
    end
    t=num2str(A(x));
    T=['frequency=' t];
    
    subplot(1,2,1)
    plot(XX,squeeze(mean(VEP_Fr(xx,:,:),2)))
    title(T);
    xlabel('Time(s)')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.YLim=[-0.1 0.1];
    ax.XLim=[0 dur_in_sec];
    
    subplot(1,2,2)
    plot(f,squeeze(mean(P_data(xx,:,:),2)),'-k')
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
    xx=xx+1;
end


clear x P ft

% check cycle time in metropsis/VEP set up

clear; clc;

% load saved raw data
load 'test_VEP_dataCycletime.mat'

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
D=0;
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
dur_in_sec=2; % the length of the stimulus presentation
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
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.YLim=[-0.1 0.1];
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


% bootstrap analysis for cycletime
[Bootstat,Bootsam]=bootstrp(100,@mean,parsed_VEP);


t=(0:dur_in_freq-1)*1/Fs;       % Time vector
f = Fs*(0:(dur_in_freq/2))/dur_in_freq; 

for x=1:size(Bootstat,1)
    ft=fft(Bootstat(x,:));
    P = abs(ft/dur_in_freq);
    P_btstrp(x,:) = P(1:dur_in_freq/2+1);
    figure(3)
    
    plot(f,P_btstrp(x,:),'-k')
    title('Bootstrapped mean data')
    xlabel('frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 500];
    %pause
    hold off
end

clear x P ft

for x=1:size(parsed_VEP,1)
    ft=fft(parsed_VEP(x,:));
    P = abs(ft/dur_in_freq);
    P_data(x,:) = P(1:dur_in_freq/2+1);
    figure(3)
    
    plot(f,P_data(x,:),'-k')
    title('Original data')
    xlabel('frequency')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[0 120];
    pause
    hold off
end

clear x P ft

parsed_VEPm=mean(parsed_VEP);
ft=fft(parsed_VEPm);
P = abs(ft/dur_in_freq);
P_dataM = P(1:dur_in_freq/2+1);
figure(3)
plot(f,P_dataM,'-k')
title('Original data mean')
xlabel('frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XLim=[0 120];


Hz16=find(f==16);
P_btstrp16=mean(P_btstrp(:,Hz16));
disp(['power 16Hz bootstrap values=' num2str(P_btstrp16)])
P_data16=mean(P_data(:,Hz16));
disp(['power 16Hz raw data averaged=' num2str(P_data16)])
disp(['power 16Hz average of raw data=' num2str(P_dataM(Hz16))])


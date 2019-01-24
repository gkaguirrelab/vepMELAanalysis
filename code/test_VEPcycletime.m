% check cycle time in metropsis/VEP set up

clear; clc;

% load saved raw data
load 'test_VEP_dataCycletime.mat'

% Bandpass filter for VEP signal
sRate=p.params.frequencyInHz;
lo=0.5; % low cut off frequency
hi=500; % high cut off frequency

d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',lo,...
    'HalfPowerFrequency2',hi,'SampleRate',sRate);
    

VEP_data=filter(d,p.response(2,:));

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
dur_in_sec=3; % the length of the stimulus presentation
dur_in_freq=dur_in_sec*sRate;
repeat=7;

for x=1:length(sync_loc)
    parsed_VEP(x,:)=VEP_data(1,sync_loc(x):sync_loc(x)+dur_in_freq);
end

clear x

XX=(1:length(parsed_VEP))/sRate;
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

% find initial peak in the mean bootstrapped data
time_in_fs=200;
for x=1:size(Bootstat,1)
    [p1,lp1]=findpeaks(Bootstat(x,1:time_in_fs),...
        sRate,'MinPeakProminence',0.025,...
        'MinPeakWidth',0.005);
    
    figure(2)
    plot(XX(:,1:time_in_fs),Bootstat(x,1:time_in_fs),'-k')
    hold on
    plot(lp1(:,1),p1(:,1),'or')
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.YLim=[-0.1 0.1];
    pause
    hold off
    
    peak1(x,:)=p1(:,1);
    loc_peak1(x,:)=lp1(:,1);
end

lp1_std=std(loc_peak1)*sRate;
disp(['standard deviation of initial peak in msec=' num2str(lp1_std)])


% function to fit sine wave to data
fit_curve=fit(XX,Bootstat(x,1:time_in_fs,));

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

% parse data by frequency
figure(2)

A=unique(trials);
VEP_Fr=zeros(length(A),repeat,length(parsed_VEP));
xx=1;
yy=1;

for x=1:length(A)
    for y=1:length(trials)
        if A(x)==trials(y)
           VEP_Fr(xx,yy,:)=parsed_VEP(y,:);
           yy=yy+1;
        end
    end
    t=num2str(A(x));
    T=['frequency=' t];
    plot(XX,squeeze(mean(VEP_Fr(xx,:,:),2)))
    title(T);
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.YLim=[-0.1 0.1];
    pause
    hold off
    
    yy=1;
    xx=xx+1;
end


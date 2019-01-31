function [p_data]=parseVEP(dur_in_sec,starttime,hi,lo)
% parseVEP
% Code to organize VEP signal and parse trials based on input, and organize
% stimuli by stimulus frequency

% VEP data is VEP.response(2,:)
% TTL pulse is VEP.response(1,:);
% timebase is VEP.timebase
clear;clc;

% Variables, change this so they are conditional
lo=0.5; % low cut off frequency for bandpass filter
hi=200; % high cut off frequency for bandpass filter
dur_in_sec=2; % the length of the VEP epuch
starttime=0.25; % if you want to start recording before (negative number) or after (positive number) the sync pulse units are in sec

% load compiled data for a single observer
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');

filenameMAT=fullfile(getpref('vepMELAanalysis','melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID],['Exp' expID '_' observerID 'compiled.mat']);

open(filenameMAT);
VEP_main=ans.VEP;

clear ans

% concatonate data across sessions
TF_trials=[];
VDS=[];
VEP=[];

for x=1:length(VEP_main)
    
    vds=VEP_main(x).VDS;
    VDS=cat(1,VDS,vds);
    
    tf_trials=VEP_main(x).mtrp.TFtrials;
    TF_trials=cat(1,TF_trials,tf_trials);
end

clear vds tf_trials

for AA=1:length(VEP_main)
    VEP=VEP_main(AA).VEP;
    VEP_data=VEP.response(2,:);
    Fs=VEP.params.frequencyInHz;
    
    % Bandpass filter for VEP signal
    d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',lo,...
        'HalfPowerFrequency2',hi,'SampleRate',Fs);
    VEP_data=filter(d,VEP_data);
    clear d
    
%      % Bandstop filter for 60 Hz noise
%     d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,...
%         'HalfPowerFrequency2',61,'SampleRate',Fs);
%     VEP_data=filter(d,VEP_data);
%     clear d

    % Find timestamp of TTL pulses
    TTL=VEP.response(1,:);
    timestamp=VEP.timebase;
    y=0;
    startFs=starttime*Fs;

    for x=2:length(TTL)
        if TTL(1,x)>4 && TTL(1,(x-1))<4
            y=y+1;
            sync_pulse(y)=timestamp(1,x+startFs);
            sync_loc(y)=x+startFs;
        end
    end

    % if the VEP recording started after the first sync pulse
    if length(sync_pulse)==35
        disp(['err sync pulse session:' num2str(AA)])
        pause
        sync_pulse=cat(2,timestamp(1,1),sync_pulse);
        sync_loc=cat(2,1,sync_loc);
    end

    clear x y

    % parse VEP data
    dur_in_freq=dur_in_sec*Fs;

    for x=1:length(sync_loc)
        parsed_vep(x,:)=VEP_data(1,sync_loc(x):sync_loc(x)+dur_in_freq);
    end

    clear x

    XX=(1:length(parsed_vep))/Fs;
    YY=mean(parsed_vep,1);
    
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
    pause
    
%     for x=1:(size(parsed_vep,1))
%         plot(XX,parsed_vep(x,:))
%         hold off
%         pause
%     end
    
    parsed_VEP(AA,:,:)=parsed_vep;
    
    clear x sync_pulse sync_loc
end

clear AA x y ax D TTL timestamp parsed_vep VEP_main repeat


% parse data by frequency

A=unique(TF_trials);
yy=1;

f = Fs*(0:(dur_in_freq/2))/dur_in_freq; 

% creates 4D matrix where VEP_Fr is session, TF, repeat, time
for w=1:size(parsed_VEP,1)
    for x=1:length(A)
        for y=1:length(TF_trials)
            if A(x)==TF_trials(w,y)
               VEP_Fr(w,x,yy,:)=parsed_VEP(w,y,:);
               VDS_Fr(w,x,yy)=VDS(w,y);
               yy=yy+1;
            end
        end
        


        figure(2)
        plot(XX,squeeze(squeeze(mean(VEP_Fr(w,x,:,:),3))))
        title(['frequency=' num2str(A(x))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 dur_in_sec];

        yy=1;
    end
end



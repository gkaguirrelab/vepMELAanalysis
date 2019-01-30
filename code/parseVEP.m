% parseVEP
% Code to organize VEP signal and parse trials based on input, and organize
% stimuli by stimulus frequency

% VEP data is VEP.response(2,:)
% TTL pulse is VEP.response(1,:);
% timebase is VEP.timebase

% Variables
lo=0.5; % low cut off frequency for bandpass filter
hi=200; % high cut off frequency for bandpass filter
dur_in_sec=1.75; % the length of the VEP epuch
starttime=0.25; % if you want to start recording before (negative number) or after (positive number) the sync pulse units are in sec

% load compiled data for a single observer
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');

filenameMAT=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_VEP/'...
         'Exp_' expID '/Subject_' observerID '/Exp' expID '_' observerID 'compiled.mat']);

open(filenameMAT);
VEP_main=ans.VEP;

clear expID observerID ans

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
    

    % Bandpass filter for VEP signal
    Fs=VEP.params.frequencyInHz;

    d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',lo,...
        'HalfPowerFrequency2',hi,'SampleRate',Fs);
    
    VEP_data=VEP.response(2,:);

    VEP_data=filter(d,VEP_data);


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
figure(2)

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
               
               % Fourier transform
                ft=fft(squeeze(parsed_VEP(w,y,:)));
                P = abs(ft/dur_in_freq);
                P_data(w,x,yy,:) = P(1:dur_in_freq/2+1);
                
                % select power spectra for each frequency
                temp=find(round(A(x))==f);
                P_dataFr(w,x,yy,:)=P_data(w,x,yy,temp);

               yy=yy+1;
            end
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

temp=find(round(A(6))==f);
Noise60=squeeze(P_data(:,:,:,temp));

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


figure(4)
% Plot mean Visual discomfort data
subplot(2,1,1)
p_dataFr=squeeze(p_dataFr);
P_Fm=mean(p_dataFr,2);
P_Fm2=mean(p_dataFr(6,:),2)-mean(noise60);
P_Fstd=std(p_dataFr,[],2);
errorbar(A,P_Fm,P_Fstd,'-ok')
hold on
errorbar(A(6),P_Fm2,P_Fstd(6),'-ob')
ylabel('power spectra for stimulus frequency')
xlabel('frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';

% Plot mean Visual discomfort data
subplot(2,1,2)
VDSm=mean(vds_Fr,2);
VDSstd=std(vds_Fr,[],2);
errorbar(A,VDSm,VDSstd,'-ok')
ylabel('visual discomfort scale')
xlabel('frequency')
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.XScale='log';
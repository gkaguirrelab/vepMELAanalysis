function [parsedVEPdata]=parseVEP(VEP_main, varargin)
% parseVEP
% Code to organize VEP signal and parse trials based on input, and organize
% stimuli by stimulus frequency

% VEP data is VEP.response(2,:)
% TTL pulse is VEP.response(1,:);
% timebase is VEP.timebase

%% Parse input
p = inputParser;
p.addParameter('lo_freq',0.5,@isnumeric); % low frequency for bandpass filter to remove noise from VEP data
p.addParameter('hi_freq',200,@isnumeric); % high frequency for bandpass filter to remove noise from VEP data
p.addParameter('dur_in_sec',2,@isnumeric); % length (in seconds) of VEP epoch to take from each trial
p.addParameter('starttime',0,@isnumeric);
p.addParameter('plot_sessions',false,@islogical);
p.addParameter('plot_all',false,@islogical);
p.addParameter('bandstop60',false,@islogical);
p.parse(varargin{:});

TF_trials=[];
VDS=[];
VEP=[];

% concatonate VDS and TF trial data across sessions
for x=1:length(VEP_main)
    
    vds=VEP_main(x).vds;
    VDS=cat(1,VDS,vds);
    
    tf_trials=VEP_main(x).mtrp.TFtrials;
    TF_trials=cat(1,TF_trials,tf_trials);
end

for AA=1:length(VEP_main)
    VEP=VEP_main(AA).vepDataStruct;
    VEP_data=VEP.response(2,:);
    Fs=VEP.params.frequencyInHz;
    
    % Bandpass filter for VEP signal
    d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',p.Results.lo_freq,...
        'HalfPowerFrequency2',p.Results.hi_freq,'SampleRate',Fs);
    VEP_data=filter(d,VEP_data);
    clear d

    if p.Results.bandstop60==1
        % Bandstop filter for 60Hz noise in VEP signal
        d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,...
            'HalfPowerFrequency2',61,'SampleRate',Fs);
        VEP_data=filter(d,VEP_data);
        clear d
    end
    
    % Find timestamp of TTL pulses
    TTL=VEP.response(1,:);
    timestamp=VEP.timebase;
    y=0;
    startFs=p.Results.starttime*Fs;

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
        sync_pulse=cat(2,timestamp(1,1),sync_pulse);
        sync_loc=cat(2,1,sync_loc);
    end

    clear x y

    % parse VEP data by trial
    dur_in_freq=p.Results.dur_in_sec*Fs;

    for x=1:length(sync_loc)
        parsed_vep(x,:)=VEP_data(1,sync_loc(x):sync_loc(x)+dur_in_freq);
    end

    clear x

    XX=(1:length(parsed_vep))/Fs;
    YY=mean(parsed_vep,1);
    
    if p.Results.plot_sessions==1
        figure(1)
        subplot(3,1,AA)
        plot(XX,YY,'-k')
        title(['Session:' num2str(VEP_main(AA).expParam.sessionID) ', average across frequencies']);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 p.Results.dur_in_sec];
    end
    
    parsed_VEP(AA,:,:)=parsed_vep;
    
    clear x sync_pulse sync_loc
end

% parse data by frequency

A=unique(TF_trials);
yy=1;

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
        
        if p.Results.plot_sessions==1
            % Plot mean frequency data for each session
            figure(w+1)
            subplot(3,2,x)
            plot(XX,squeeze(squeeze(mean(VEP_Fr(w,x,:,:),3))))
            title(['Session:' num2str(VEP_main(w).expParam.sessionID) ', frequency=' num2str(A(x))]);
            xlabel('Time(s)')
            ax=gca;
            ax.TickDir='out';
            ax.Box='off';
            ax.YLim=[-0.1 0.1];
            ax.XLim=[0 p.Results.dur_in_sec];
        end

        yy=1;
    end
end

% Concatonate data across sessions
vep_Fr=[];
vds_Fr=[];
for ZZ=1:size(VEP_Fr,1)
    vep_Fr=cat(2,vep_Fr,squeeze(VEP_Fr(ZZ,:,:,:)));
    vds_Fr=cat(2,vds_Fr,squeeze(VDS_Fr(ZZ,:,:)));
end

if p.Results.plot_all==1
    % Plot mean frequency data across sessions
    figure(7)
    for x=1:length(A)
        subplot(3,2,x)
        plot(XX,squeeze(mean(vep_Fr(x,:,:),2)))
        title(['frequency=' num2str(A(x))]);
        xlabel('Time(s)')
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.YLim=[-0.1 0.1];
        ax.XLim=[0 p.Results.dur_in_sec];
       
    end
        % Plot mean Visual discomfort data
        figure(8)
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

end

parsedVEPdata.parsed_VEP=parsed_VEP;
parsedVEPdata.VEP_Fr=VEP_Fr;
parsedVEPdata.vep_Fr=vep_Fr;
parsedVEPdata.VDS_Fr=VDS_Fr;
parsedVEPdata.vds_Fr=vds_Fr;
parsedVEPdata.dur_in_freq=dur_in_freq;
parsedVEPdata.Fs=Fs;

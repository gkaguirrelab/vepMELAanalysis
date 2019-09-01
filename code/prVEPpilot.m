% variables
nTrials=60; % number of trials in file
% calls the experiment based on users input
expID='PVEP1';
observerID='CPG';
sessionID=[1,2]; 

% Establish directory for saving. The path is defined by experiment ID.
savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID]);
if ~exist(savePath,'dir')
    mkdir(savePath);
end

filenameComp=fullfile(savePath,['Exp' expID '_' observerID 'prvep.mat']);
     

    for x=1:length(sessionID)
        % Get VEP data file
         filenameMAT=fullfile(getpref('vepMELAanalysis', 'melaDataPath'), 'Experiments','OLApproach_VEP',...
             ['Exp_' expID], ['Subject_' observerID], ['Exp' expID '_' observerID num2str(sessionID(x)) '.mat']);

          vep=open(filenameMAT);

      % Get metropsis data, load metropsis JSON file
          % Load and read JSON file
    filenameJSON=fullfile(getpref('vepMELAanalysis', 'mtrpDataPath'),['Exp_' expID],...
        ['Subject_' observerID], [observerID '_' num2str(sessionID(x)) '.lpkmx']);  
        [vep.mtrp]=parse_mtrp_json(filenameJSON,nTrials);
        
        VEP(x)=vep;
    end
    
    clear vep trial_dur nTrials x
    
[parsedVEPdata]=parse_vep(VEP);

vep_dur=(0:length(parsedVEPdata.VEP_dataC))/parsedVEPdata.Fs;
vep_dur=vep_dur(1:end-1);
vepM=squeeze(nanmean(parsedVEPdata.VEP_dataC,2));

figure(1)
hold on
plot(vep_dur,vepM(1,:),'Color',[0.8 0.8 0.8])
plot(vep_dur,vepM(2,:),'Color',[0.5 0.5 0.5])
plot(vep_dur,vepM(3,:),'k')
ax=gca;ax.Box='off';ax.TickDir='out';
    
for x=1:size(parsedVEPdata.VEP_dataC,1)
    for y=1:size(parsedVEPdata.VEP_dataC,2)
        vep=squeeze(squeeze(parsedVEPdata.VEP_dataC(x,y,:)))';
%         vep=vep-vepM(x,:); % induced VEP
%         [psd(x,y,:),freqs]=pmtm(vep',3,length(vep),parsedVEPdata.Fs);
        [psd(x,y,:),freqs]=pwelch(vep',[],[],[],parsedVEPdata.Fs);
    end
end

clear vep

PSD=squeeze(nanmean(psd,2));

figure(2)
hold on
plot(freqs',PSD(1,:),'Color',[0.8 0.8 0.8])
plot(freqs',PSD(2,:),'Color',[0.5 0.5 0.5])
plot(freqs',PSD(3,:),'k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 150];ax.YScale='log';


for x=1:size(parsedVEPdata.VEP_dataC,1)
    for y=1:size(parsedVEPdata.VEP_dataC,2)
        vep=squeeze(parsedVEPdata.VEP_dataC(x,y,:))';
%         vep=vep-vepM(x,:); % induced VEP
        [S(x,y,:,:),F,T]=spectrogram(vep,[],[],[],parsedVEPdata.Fs);
    end
end

S=squeeze(nanmean(S,2));

figure(3)
subplot(1,3,1)
helperCWTTimeFreqPlot(squeeze(S(1,:,:)),T,F,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];

subplot(1,3,2)
helperCWTTimeFreqPlot(squeeze(S(2,:,:)),T,F,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];

subplot(1,3,3)
helperCWTTimeFreqPlot(squeeze(S(3,:,:)),T,F,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];


% pattern reversal VEPs
[prVEP]=parse_prVEP(parsedVEPdata.VEP_dataC,parsedVEPdata.Fs);
prVEPm=squeeze(nanmean(prVEP,2));
PRvep_dur=(0:length(prVEP))/parsedVEPdata.Fs;
PRvep_dur=PRvep_dur(1:end-1);

for x=1:size(prVEPm,1)
    temp=nanmean(prVEPm(x,1:100));
    prVEPm(x,:)=prVEPm(x,:)-temp*ones(size(prVEPm(x,:))); 
end

figure(4)
hold on
plot(PRvep_dur,prVEPm(1,:),'Color',[0.8 0.8 0.8])
plot(PRvep_dur,prVEPm(2,:),'Color',[0.5 0.5 0.5])
plot(PRvep_dur,prVEPm(3,:),'k')
ax=gca;ax.Box='off';ax.TickDir='out';
    
for x=1:size(prVEP,1)
    for y=1:size(prVEP,2)
        vep=squeeze(squeeze(prVEP(x,y,:)))';
%         vep=vep-prVEPm(x,:); % induced VEP
        [PRpsd(x,y,:),freqs]=pwelch(vep',[],[],[],parsedVEPdata.Fs);
%         [PRpsd(x,y,:),freqs]=pmtm(vep,3,length(vep),parsedVEPdata.Fs);
    end
end

clear vep

prPSD=squeeze(nanmean(PRpsd,2));

figure(5)
hold on
plot(freqs',prPSD(1,:),'Color',[0.8 0.8 0.8])
plot(freqs',prPSD(2,:),'Color',[0.5 0.5 0.5])
plot(freqs',prPSD(3,:),'k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 150];ax.YScale='log';


for x=1:size(prVEP,1)
    for y=1:size(prVEP,2)
        vep=squeeze(prVEP(x,y,:))';
%         vep=vep-prVEPm(x,:); % induced VEP
        [prS(x,y,:,:),prF,prT]=spectrogram(vep,[],[],[],parsedVEPdata.Fs);
    end
end

prS=squeeze(nanmean(prS,2));

figure(6)
subplot(1,3,1)
helperCWTTimeFreqPlot(squeeze(prS(1,:,:)),prT,prF,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];

subplot(1,3,2)
helperCWTTimeFreqPlot(squeeze(prS(2,:,:)),prT,prF,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];

subplot(1,3,3)
helperCWTTimeFreqPlot(squeeze(prS(3,:,:)),prT,prF,'surf','STFT for VEP signal','seconds','Hz')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 80];

%% local functions
function [mtrp]=parse_mtrp_json(filenameJSON,nTrials)

    json=fileread(filenameJSON);

    % Use json file to find the start time, observer ID, group, and session number of the experiment
    t=strfind(json,'Starting Time = ')+length('Starting Time = ');
    s=strfind(json,'Session #')+length('Session #');
    o=strfind(json,'Subject: ')+length('Subject: ');
    g=strfind(json,'Group: ')+length('Group: ');

    mtrp.date=json(t:t+9);
    mtrp.starttime=json(t+11:t+18);
    mtrp.group=json(g:g+3);
    mtrp.observer=json(o:o+2);

    if isempty(str2num(json(s+1)))
        mtrp.session=num2str(json(s));
    else
        mtrp.session=num2str(json(s:s+1));
    end


    % find the contrast values

    k=strfind(json,'contrast]" = ')+length('contrast]" = '); % the number of characters that the numeric value of contrast is after the searched string
    k=k(end-(nTrials-1):end);
    for x=1:length(json)
        for y=1:length(k)
            if x==k(y) && ~isempty(str2num(json(k(y))))
                trial_C(y)=str2num(json(k(y):k(y)+1));
            else if x==k(y) && isempty(str2num(json(k(y))))
                   trial_C(y)=str2num(json(k(y)+1:k(y)+3)); 
                end
            end
        end
    end
    
    for a=1:length(trial_C)
        if trial_C(a)==10
            trial_C(a)=100;
        end
    end
    
    mtrp.Ctrials=trial_C;
end




function [parsedVEPdata]=parse_vep(VEP)

lo_freq=0.5;
hi_freq=200;
dur_in_sec=1;
C_trials=[];
VEP_data=[];

for x=1:length(VEP)
    
    % concatonate contrast values for trials
    temp=VEP(x).mtrp.Ctrials';
    C_trials=cat(1,C_trials,temp);
    
    % filter and concatonate parsed VEP data for each trial
    Fs=VEP(x).vepDataStruct.params.frequencyInHz;
    ttl=VEP(x).vepDataStruct.response(1,:);
    timestamp=VEP(x).vepDataStruct.timebase;
    y=0;
    % Find timestamp of TTL pulses
    for xx=2:length(ttl)
        if ttl(1,xx)>4 && ttl(1,xx)<6 && ttl(1,(xx-1))<4
            y=y+1;
            sync_pulse(y)=timestamp(1,xx);
            sync_loc(y)=xx;
        end
    end
    
    vep=VEP(x).vepDataStruct.response(2,:);
     % Bandpass filter for VEP signal
    d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',lo_freq,...
        'HalfPowerFrequency2',hi_freq,'SampleRate',Fs);
    vep=filter(d,vep);
    clear d
    % Bandstop filter for 60Hz noise in VEP signal
    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,...
        'HalfPowerFrequency2',61,'SampleRate',Fs);
    vep=filter(d,vep);
    clear d
    % Bandstop filter for 120Hz noise in VEP signal
    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',119,...
        'HalfPowerFrequency2',121,'SampleRate',Fs);
    vep=filter(d,vep);
    clear d
    
    % parse VEP data by trial
    dur_in_freq=dur_in_sec*Fs;

    for xxx=1:length(sync_pulse)
        parsedVEPsession(xxx,:)=vep(1,sync_loc(xxx):sync_loc(xxx)+dur_in_freq);
    end
    
    VEP_data=cat(1,VEP_data,parsedVEPsession);
end

clear temp*

% parse VEP data by contrast
contrast_value=[0 50 100];
for X=1:length(contrast_value)
    temp=find(C_trials==contrast_value(X));
    VEP_dataC(X,:,:)=VEP_data(temp,:);
end

parsedVEPdata.VEP_dataC=VEP_dataC;
parsedVEPdata.Fs=Fs;

end




function [prVEP]=parse_prVEP(vep_data,Fs)

    prAll=[];
    for x=1:size(vep_data,1)
        for y=1:size(vep_data,2)
            temp=squeeze(vep_data(x,y,:))';
            start1=0.25*Fs;
            start2=0.5*Fs;
            start3=0.75*Fs;
            dur=0.25*Fs;
            pr1=temp(:,start1:start1+dur);
            pr2=temp(:,start2:start2+dur);
            pr3=temp(:,start3:start3+dur);
            prAll=cat(1,prAll,pr1,pr2,pr3);
            clear pr1 pr2 pr3  
        end
        prVEP(x,:,:)=prAll;
        prAll=[];
    end
    
end

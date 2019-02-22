function [VEP]=compileVEPdata()
% A function that compiles data for a single subject across sessions, and
% organizes experiment parameters, audio data, and VEP data
%
% Syntax:
%  [metVEPdataStruct]=compileVEPdata()
%
% Description:
%	This function organizes experimental data from the VEP/metropsis experiments for a
%	single subject across all the sessions they ran.
%
%
% Output: saves structure VEP
%   VEP                   - A structure that contains VEP data, TTL pulse,
%                           and timebase for both analog signals
%   audioRec              - A structure that contains mic data, and 
%                           sampling rate (Fs)
%   expParam              - contains observer ID, experiment ID, session 
%                           ID, data that was recorded from the VEP computer
%   mtrp                  - structure that contains all the metropsis data. 
%                           This includes observer ID, group (MWVA or HA free), 
%                           session, and temporal frequency for the stimuli 
%                           presented (TFtrials)
%   VDS                   - Visual discomfort scale values for the 36 trials

% variables
nTrials=25; % number of trials in file
trial_dur=5; % trial duration in seconds
% calls the experiment based on users input
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');
sessionID=input('session ID:'); 

% Establish directory for saving. The path is defined by experiment ID.
savePath = fullfile(getpref('vepMELAanalysis', 'melaAnalysisPath'),'experiments',...
    'vepMELAanalysis',['Exp_' expID]);
if ~exist(savePath,'dir')
    mkdir(savePath);
end

filenameComp=fullfile(savePath,['Exp' expID '_' observerID 'compiled.mat']);
     
Z=[];
if exist(filenameComp)==2
         disp('This dataset already exists!')
         Z=input('Do you want to continue? (Y/N):');
else
    Z=('Y');
end

if Z==('N') && exist(filenameComp)==2
    disp('compileVEPdata aborted')
else
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

        % Get audio data
        [vep.vds]=parseAudioData(vep,nTrials,trial_dur);
        
        VEP(x)=vep;
    end
    
       save(filenameComp,'VEP')
if exist(filenameComp)==2
     disp('compiled data saved')
     pause
    clear; clc;
else
     disp('error, file did not save')
end 
    
end
end

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


    % find the temporal frequencies

    k=strfind(json,'tf]" = ')+length('tf]" = '); % the number of characters that the numeric value of TF is after the searched string
    k=k(end-(nTrials-1):end);
    for x=1:length(json)
        for y=1:length(k)
            if x==k(y) && ~isempty(str2num(json(k(y))))
                trial_TF(y)=str2num(json(k(y):k(y)+1));
            else if x==k(y) && isempty(str2num(json(k(y))))
                   trial_TF(y)=str2num(json(k(y)+1:k(y)+3)); 
                end
            end
        end
    end
    
    for a=1:length(trial_TF)
        if trial_TF(a)==3.2
            trial_TF(a)=3.25;
        end

        if trial_TF(a)==1.6
            trial_TF(a)=1.625;
        end
    end
    
    mtrp.TFtrials=trial_TF;

end



function [VDS]=parseAudioData(VEP,nTrials,trial_dur)
    audioRec=VEP.audioRec;

    % Check audio files
    X=((1:1:length(audioRec.data))./audioRec.Fs);
    parse_audio=X(1:trial_dur*audioRec.Fs:end);
    figure(1)
    plot(X,audioRec.data)
    hold on
    plot(parse_audio,0.05*ones(1,length(parse_audio)),'or')
    hold off
    shift_time=input('shift the time=');
    parse_audio=X(audioRec.Fs*shift_time:trial_dur*audioRec.Fs:end);

    for tt = 1:nTrials
        firstTrialIndex = parse_audio(tt)*audioRec.Fs;
        secondTrialIndex = parse_audio(tt+1)*audioRec.Fs;
        trialAudio = audioRec.data(firstTrialIndex:secondTrialIndex);
        [ firstTimePoint, secondTimePoint ] = grabRelevantAudioIndicesVEP(trialAudio, audioRec.Fs);

        % index that starts the numerical rating
        firstIndex = firstTimePoint*audioRec.Fs;

        % index that ends the numerical rating
        secondIndex = secondTimePoint*audioRec.Fs;

        if firstIndex<1
            firstIndex=1;
        end

        if isempty(firstIndex) || isempty(secondIndex) || secondIndex>length(trialAudio)
            sound(trialAudio, audioRec.Fs);
        else
            sound(trialAudio(firstIndex:secondIndex), audioRec.Fs);
        end

        vds=input('visual discomfort scale:');

        if isnan(vds)
            sound(trialAudio, audioRec.Fs);
            vds=input('visual discomfort scale:');
        end

        VDS(tt)=vds;
    end
end
% Main script that runs metropsis VEP analysis

% calls the experiment based on users input
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');
sessionID=input('session ID:');

filenameComp=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_VEP/'...
         'Exp_' expID '/Subject_' observerID '/Exp' expID '_' observerID 'compiled.mat']);
     
Z=[];
if exist(filenameComp)==2
         disp('This dataset already exists!')
         Z=input('Do you want to continue? (Y/N):');
end

if Z==('N') && exist(filenameComp)==2
    disp('runVEPanalysis aborted')
else
    
    for x=1:length(sessionID)

        % Get VEP data
         filenameMAT=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_VEP/'...
             'Exp_' expID '/Subject_' observerID '/Exp' expID '_' observerID num2str(sessionID(x)) '.mat']);

          vep=open(filenameMAT);

      % Get MTRP data
        [vep.mtrp]=parse_mtrp_json(expID,observerID,num2str(sessionID(x)));
        VEP(x)=vep;

        % Check timestamp for metropsis file and VEP data - I need to save time
        % for VEP data
    end

    %% Get audio files, parse them, and record their numeric value
    nTrials=36; % number of trials in file
    trial_dur=6; % trial duration in seconds
    
    for xx=1:length(VEP)
        audioRec=VEP(xx).audioRec;
        
        % Check audio files
            X=((1:1:length(audioRec.data))./audioRec.Fs);
            parse_audio=X(1:trial_dur*audioRec.Fs:end);
            figure(1)
            plot(X,audioRec.data)
            hold on
            plot(parse_audio,0.05*ones(1,length(parse_audio)),'or')
            hold off
            pause
            
            
            
        for tt = 1:nTrials
            firstTrialIndex = parse_audio(tt)*audioRec.Fs;
            secondTrialIndex = parse_audio(tt+1)*audioRec.Fs;
            trialAudio = audioRec.data(firstTrialIndex:secondTrialIndex);
            [ firstTimePoint, secondTimePoint ] = grabRelevantAudioIndices(trialAudio, audioRec.Fs);

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
        VEP(xx).VDS=VDS;
        clear VDS parse_audio
    end

% save compiled data


     save(filenameComp,'VEP')
     if exist(filenameComp)==2
         disp('compiled data saved')
         pause
        clear; clc;
     else
         disp('error, file did not save')
     end
end     
%% related functions
function [mtrp]=parse_mtrp_json(expID,observerID,sessionID)
    
    % Load and read JSON file
    filenameJSON=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MTRP_data/Exp_' expID '/Subject_' observerID '/' observerID '_' sessionID '.lpkmx']); 
    json=fileread(filenameJSON);


    % Use json file to find the start time, observer ID, group, and session number of the experiment
    t=strfind(json,'Starting Time =')+16;
    s=strfind(json,'Session #')+9;
    o=strfind(json,'Subject:')+9;
    g=strfind(json,'Group:')+7;

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

    k=strfind(json,'tf]" =')+7; % 8 is the number of characters that the numeric value of TF is after the searched string
    k=k(end-35:end);
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






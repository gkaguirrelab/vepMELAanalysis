% Main script that runs metropsis VEP analysis
clear;clc;

% calls the experiment based on users input
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');
sessionID=input('session ID:','s');

%% get metropsis data

[mtrp]=parse_mtrp_xml(expID,observerID,sessionID);
%% Get VEP data
 filenameMAT=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_VEP/'...
     'Exp_' expID '/Subject_' observerID '/Exp' expID '_' observerID sessionID '.mat']);
 
 VEP=open(filenameMAT);

 clear filenameMAT expID observerID sessionID

%% Compare metropsis experiment start time to VEP start time to ensure the files are matched up properly



%% related functions
function [mtrp]=parse_mtrp_xml(expID,observerID,sessionID)
    
    % Load and read XML file
    filenameXML=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MTRP_data/Exp_' expID '/Subject_' observerID '/' observerID '_' sessionID '.dpkmx']); 
    xml=fileread(filenameXML);
    
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
  
    
    % find the temporal frequencies for each trial using XML file
 
    
    k=strfind(xml,'[tf]</key>')+25; % 25 is the number of characters that the numeric value of TF is after the searched string
    for x=1:length(xml)
        for y=1:length(k)
            if x==k(y) && ~isempty(str2num(xml(k(y)+1)))
                trial_TF(y)=str2num(xml(k(y):k(y)+1));
            else if x==k(y) && isempty(str2num(xml(k(y)+1)))
                   trial_TF(y)=str2num(xml(k(y))); 
                end
            end
        end
    end
    
    mtrp.TFtrials=trial_TF;
    
    % find channel using XML file
    S=strfind(xml,'[Scone]</key>')+24;
    
       
    switch str2num(xml(S))
        case -1
            mrtp.channel=char('(L+M)-S');
        case 0
            mrtp.channel=char('L-M');
        case 1
            mtrp.channel=char('Lum');
    end
    
end






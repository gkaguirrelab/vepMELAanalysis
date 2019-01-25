% Main script that runs metropsis VEP analysis
clear;clc;

% Get VEP data
expID=input('experiment ID:','s');
observerID=input('observer ID:','s');
sessionID=input('session ID:','s');

 filenameXML=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_VEP/'...
     expID '/Subject_' observerID '/Exp' expID '_' observerID sessionID '.mat']);

% get metropsis data
%[mtrp]=parse_mtrp_xml(expID,observerID,sessionID);



%% related functions
function [mtrp]=parse_mtrp_xml(expID,observerID,sessionID)

% Read and decode metropsis data witht he suffix .dpkmx, which is in an xml file

    filenameXML=fullfile(['/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MTRP_data/Exp_' expID '/Subject_' observerID '/' observerID '_' sessionID '.dpkmx']); 
    xml=fileread(filenameXML);
    k=strfind(xml,'[tf]</key>')+25; % 25 is the number of characters that the numeric value of TF is after the searched string

    for x=1:length(xml)
        for y=1:length(k)
            if x==k(y)
                trial_TF(y)=str2num(xml(k(y)));
            end
        end
    end
    
    k=strfind(xml,'[Lcone]</key>');
    
    mtrp.temporalF=trial_TF;
    
end






function [ttf_fit, TemporalFrequency_fit,params]=getTTFfits_VDS(flicker_response,stimulusFreqHz,x0)
    % Fit the Watson model to data
   
    % Find the maximum interpolated VEP response
    stimulusFreqHzFine=logspace(0,log10(max(stimulusFreqHz)),100);
    
    % Scale the x vector so that the max is 1
    scaledResponse=flicker_response./max(flicker_response);
%     x0=cat(2,x0,max(scaledVEP));
    myObj=@(p)sqrt(sum((scaledResponse-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    x0_lb=x0-[1 0.5 0.5];
    x0_ub=x0+[1 0.5 0.5];
    params=fmincon(myObj,x0,[],[],[],[],x0_lb,x0_ub);
    
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine,params).*max(flicker_response);
    TemporalFrequency_fit=stimulusFreqHzFine;
end
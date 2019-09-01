function [ttf_fit, TemporalFrequency_fit,params]=getTTFfits(VEPresponse,stimulusFreqHz,x0)
    % Fit the Watson model to data
    % Adjust the VEP response to deal with negative values
    minVEP=min(VEPresponse);
    if minVEP<0
        scaledVEP=VEPresponse-minVEP;
    else
        scaledVEP=VEPresponse;
        minVEP=0;
    end
    % Find the maximum interpolated VEP response
    stimulusFreqHzFine=logspace(0,log10(max(stimulusFreqHz)),100);
    splineInterpolatedMax=max(spline(stimulusFreqHz,VEPresponse,stimulusFreqHzFine));
    
    % Scale the x vector so that the max is 1
    scaledVEP=scaledVEP./max(scaledVEP);
%     x0=cat(2,x0,max(scaledVEP));
    myObj=@(p)sqrt(sum((scaledVEP-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    x0_lb=x0-[1 0.5 0.5];
    x0_ub=x0+[1 0.5 0.5];
    params=fmincon(myObj,x0,[],[],[],[],x0_lb,x0_ub);
    
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine,params).*splineInterpolatedMax+minVEP;
    TemporalFrequency_fit=stimulusFreqHzFine;
end
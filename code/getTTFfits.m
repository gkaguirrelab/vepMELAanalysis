function [ttf_fit, TemporalFrequency_fit,params]=getTTFfits(VEPresponse,stimulusFreqHz,x0,varargin)
    % Fit the Watson model to data
    q = inputParser;
    q.addParameter('x0_lb',x0-[1 0.5 0.5],@isnumeric);
    q.addParameter('x0_ub',x0+[1 0.5 0.5],@isnumeric);
    
    q.parse(varargin{:});
    
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
    myObj=@(p)sqrt(sum((scaledVEP-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    params=fmincon(myObj,x0,[],[],[],[],q.Results.x0_lb,q.Results.x0_ub);
    
    stimulusFreqHzFine2=logspace(0,log10(max(stimulusFreqHz)+20),100);
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine2,params).*max(VEPresponse)+minVEP;
    TemporalFrequency_fit=stimulusFreqHzFine2;
end
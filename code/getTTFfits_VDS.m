function [ttf_fit, TemporalFrequency_fit,params]=getTTFfits_VDS(flicker_response,stimulusFreqHz,x0,varargin)
    % Fit the Watson model to data

    q = inputParser;
    q.addParameter('x0_lb',x0-[1 0.5 0.5],@isnumeric);
    q.addParameter('x0_ub',x0+[1 0.5 0.5],@isnumeric);

    q.parse(varargin{:});

    % Find the maximum interpolated VEP response
    stimulusFreqHzFine=logspace(0,log10(max(stimulusFreqHz)),100);
    
    % Scale the x vector so that the max is 1
    scaledResponse=flicker_response./max(flicker_response);
    myObj=@(p)sqrt(sum((scaledResponse-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    params=fmincon(myObj,x0,[],[],[],[],q.Results.x0_lb,q.Results.x0_ub);
    
    ttf_fit=watsonTemporalModelvep(stimulusFreqHzFine,params).*max(flicker_response);
    TemporalFrequency_fit=stimulusFreqHzFine;
end
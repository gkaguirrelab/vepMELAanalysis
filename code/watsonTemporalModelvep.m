function y = watsonTemporalModelvep(frequenciesToModel, params, params_centerAmplitudeIn, filterOrder)
% Beau Watson's 1986 center-surround neural temporal sensitivity model
%
% Syntax:
%  y = watsonTemporalModelvep(frequenciesHz, params, params_centerAmplitude)
%
% Description:
%	Calculates the two-component (center-surround) Watson temporal model
%   The parameters (p) defines both the "center" and the "surround" linear
%   filter components. The entire model is the difference between these two
%   filter components.
%
%	The model is expressed in Eq 45 of:
%
%       Watson, A.B. (1986). Temporal sensitivity. In Handbook of Perception
%       and Human Performance, Volume 1, K. Boff, L. Kaufman and
%       J. Thomas, eds. (New York: Wiley), pp. 6-1-6-43.
%
%   Note that there is a typo in original manuscript. Equation 45 should
%   be:
%
%       H(frequenciesHz) = a[H1(frequenciesHz) - bH2(frequenciesHz)]
%
%   where a and b are scale factors. We have modified the implementation of
%   scale factors here.
%
%   Additionally, Watson (1986) gives the time-constant of the model in
%   units of milliseconds, but we find that, to reproduce the presented
%   figures, this time-constant is converted at some point to units of
%   seconds prior to its entry into the equations.
%
%   The model is the difference between two linear impulse response
%   filters, each of which is themselves a cascade of low-pass filters. The
%   number of filters in the cascade (the "filterOrder") is set
%   empirically. Center and surround orders of "9" and "10" are presented
%   in (e.g.) Figure 6.5 of Watson (1986).
%
%   The centerAmplitude parameter input to the model is optional. If
%   provided, the model will use this value in the calculation. If not
%   provided, the routine will determine a value for the centerAmplitude
%   that causes the output of the model to have a maximum value of unity
%   across a hard-coded range of plausible frequency inputs.
%
%   Note that the model can only return positive amplitudes of response,
%   although some empirical data may have a negative amplitude (e.g., BOLD
%   fMRI data of a low temporal frequency stimulus as compared to a blank
%   screen). To handle this circumstance, it is recommended that data be
%   offset so that the minimum amplitude value is zero prior to fitting.
%
% Inputs:
%   frequenciesToModel    - 1xn vector that provides the stimulus
%                           frequencies in Hz for which the model will be
%                           evaluated
%   params                - 1x3 vector of model parameters:
%                             tau - time constant of the center filter (in
%                                   msecs)
%                           kappa - multiplier surround time-constant
%                            zeta - multiplier of the surround amplitude
%                                   the surround filter
%   params_centerAmplitude - Scalar. Optional. The fourth parameter of the
%                           model.
%   filterOrder           - 1x2 vector. Optional provides the order of the
%                           center and surround filters respetively. If not
%                           provided, defaults to [9 10]
%
% Outputs:
%   y                     - 1xn vector of modeled amplitude values.
%
% Examples:
%{
    % Demonstrate basic output of the model
    freqHz = logspace(0,log10(64),100);
    params = [2 2 1];
    y = watsonTemporalModelvep(freqHz,params);
    semilogx(freqHz,y,'-k');    
%}
%{
    % Fit the Watson model to some empirical data
    stimulusFreqHz = [0.5 1 2 4 8 16 32 64];
    pctBOLDresponse = [-0.1050   -0.0380    0.1270    0.1800    0.2970    0.5110    0.4700   -0.1020];
    % Adjust the BOLD response to deal with negative values
    minBOLD = min(pctBOLDresponse)
    if minBOLD < 0
        scaledBOLDresponse = pctBOLDresponse - minBOLD;
    else
        scaledBOLDresponse = pctBOLDresponse;
        minBOLD = 0;
    end
    % Find the maximum interpolated BOLD response
    stimulusFreqHzFine = logspace(0,log10(64),100);
    splineInterpolatedMax = max(spline(stimulusFreqHz,scaledBOLDresponse,stimulusFreqHzFine));
    % Scale the x vector so that the max is zero
    scaledBOLDresponse = scaledBOLDresponse ./ splineInterpolatedMax;
    myObj = @(p) sqrt(sum((scaledBOLDresponse-watsonTemporalModelvep(stimulusFreqHz,p)).^2));
    x0 = [4 2 1];
    params = fmincon(myObj,x0,[],[]);
    semilogx(stimulusFreqHzFine,watsonTemporalModelvep(stimulusFreqHzFine,params).*splineInterpolatedMax+minBOLD,'-k');
    hold on
    semilogx(stimulusFreqHz, pctBOLDresponse, '*r');
    hold off
%}

% Handle optional parameters
if nargin >=3
    params_centerAmplitude = params_centerAmplitudeIn;
else
    params_centerAmplitude = [];
end

if nargin >=4
    centerFilterOrder = filterOrder(1); % Order of the center (usually fast) filter
    surroundFilterOrder = filterOrder(2); % Order of the surround (usually slow) filter
else
    % Fixed parameters (taken from Figure 6.4 and 6.5 of Watson 1986)
    centerFilterOrder = 9; % Order of the center (usually fast) filter
    surroundFilterOrder = 10; % Order of the surround (usually slow) filter
end


% Define a frequency domain in Hz over which the model is defined. The
% maximum and minimum value of the y response should be contained within
% this range for any plausible set of parameters. We hard code a log-spaced
% range of 0.1 - 200 Hz here.
freqDomain = logspace(-1,log10(200),100);

% Sanity check the frequency input
if max(frequenciesToModel)>max(freqDomain) || min(frequenciesToModel)<min(freqDomain)
    error('The passed frequency is out of range for the model');
end

% Un-pack the passed parameters
params_tau = params(1)./1000;   % Convert from msecs to secs
params_kappa = params(2);
params_zeta = params(3);

% Calculte the response. If a centerAmplitude was passed, perform the
% computation
if ~isempty(params_centerAmplitude)
    H1 = nStageLowPassFilter(params_tau,frequenciesToModel,centerFilterOrder);
    H2 = nStageLowPassFilter(params_kappa*params_tau,frequenciesToModel,surroundFilterOrder);
    y = (params_centerAmplitudeIn * H1) - (params_zeta*params_centerAmplitudeIn*H2);
else
    % Search to find the center amplitude that provides a maximum response
    % of unity, then perform the computation. This is a recursive call to
    % this function.
    myObj = @(x) abs(max(watsonTemporalModelvep(freqDomain, params, x, [centerFilterOrder surroundFilterOrder]))-1);
    params_centerAmplitudeIn = fminsearch(myObj,1);
    H1 = nStageLowPassFilter(params_tau,frequenciesToModel,centerFilterOrder);
    H2 = nStageLowPassFilter(params_kappa*params_tau,frequenciesToModel,surroundFilterOrder);
    y = (params_centerAmplitudeIn * H1) - (params_zeta*params_centerAmplitudeIn*H2);
end

% The model calculates a vector of complex values that define the Fourier
% transform of the system output. As we are implementing a model of just
% the amplitude component of a temporal transfer function, the absolute
% value of the model output is returned.
y = abs(y);

end % main function


function Hsub = nStageLowPassFilter(tau,frequenciesHz,filterOrder)
% This function implements the system response of the linear filter
% for temporal sensitivity of neural systems in Eq 42 of Watson (1986).
%
% The implemented function is the "system respone" (Fourier transform) of
% the impulse response of an nth-order filter which is of the form:
%
% h(t) = u(t) * (1/(tau*(n-1)!)) * (t/tau)^(n-1) * exp(-t/tau)
%
% tau -- Time constant of the filter (in seconds)
% frequenciesHz -- a vector of frequencies at which to realize the model
% filterOrder -- the number of low-pass filters which are cascaded in
%                the model
Hsub = (1i*2*pi*frequenciesHz*tau + 1) .^ (-filterOrder);

end
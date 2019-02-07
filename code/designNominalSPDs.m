function status = designNominalSPDs(varargin)
% Nominal primaries and SPDs for isolating post-receptoral mechanisms
%
% Syntax:
%	designNominalSPDs
%
% Description:
%   Given a calibration file from a monitor, this code will provide the RGB
%   primary settings that should serve to produce isolated stimulation of
%   the cone (L, M, S) and post-receptoral channels (LMS, L-M, S).
%
% Inputs:
%	None
%
% Outputs:
%	primaries             - Struct. The primaries and SPDs.
%
% Optional key/value pairs:
%  'calFilePath'          - Char. Full path to the .mat cal file.
%  'plotDir'              - Char. Full path to the directory in which the
%                           diagnostic plots will be saved. The directory
%                           will be created if it does not exist. If set to
%                           empty, no files will be saved.
%  'primaryHeadRoom'      - Scalar. We can enforce a constraint that we
%                           don't go right to the edge of the gamut.  The
%                           head room parameter is defined in the [0-1]
%                           device primary space.  Using a little head room
%                           keeps us a bit away from the hard edge of the
%                           device.
%  'observerAgeInYears'   - Scalar
%  'fieldSizeDegrees'     - Scalar
%  'pupilDiameterMm'      - Scalar
%

%% Clear and close
close all;

%% Parse input
p = inputParser;
p.addParameter('calFilePath',getpref('vepMELAanalysis','calFilePath'),@ischar);
p.addParameter('plotDir','~/Desktop/nominalSPDs',@ischar);
p.addParameter('validationFractionTolerance',0.001,@isscalar)
p.addParameter('primaryHeadRoom',0,@isscalar)
p.addParameter('observerAgeInYears',25,@isscalar)
p.addParameter('fieldSizeDegrees',30,@isscalar)
p.addParameter('pupilDiameterMm',2,@isscalar)
p.parse(varargin{:});



whichModel = 'human';
whichPrimaries = 'monitor';


% Obtain the cal file and SPDs of the primaries
[calPath, calFileName, ~] = fileparts(p.Results.calFilePath);
cal = LoadCalFile(calFileName,[],calPath);
S = [380 2 201];
B_primary = cal.processedData.P_device;
ambientSpd = cal.processedData.P_ambient;

% Set background to the monitor midpoint
backgroundPrimary = [0.5 0.5 0.5]';



%% Get sensitivities and set other relvant parameters
% The routines that do these computations are in the
% ContrastSplatter directory of the SilentSubstitutionToolbox. They
% provide pre-defined receptor types and compute spectral
% sensitivities using the routines provided in the Psychtoolbox.
% The routines here, however, also allow computation of fraction
% cone bleached, which may be used to adjust pigment peak optical
% density.  They can also compute photopigment variants corrected
% for filtering by blood vessels.

% Prompt user for key parameters that affect the spectral
% sensitivities.
%
% Note that we don't typically vary or pass the blood vessel
% parameters but rather simply accept the defaults used by
% GetHumanPhotoreceptorSS.  It's mainly for fun that we show
% how to do this here.


% Define photoreceptor classes that we'll consider.
% ReceptorIsolate has a few more built-ins than these.
photoreceptorClasses = {'LConeTabulatedAbsorbance', 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance', 'Melanopsin', 'Rods', ...
    'LConeTabulatedAbsorbancePenumbral', 'MConeTabulatedAbsorbancePenumbral', 'SConeTabulatedAbsorbancePenumbral'};

% Make sensitivities.  The wrapper routine is
% GetHumanPhotoreceptorSS, which is in the ContrastSplatter
% directory.  Each row of the matrix T_receptors provides the
% spectral sensitivity of the photoreceptor class in the
% corresponding entry of the cell array photoreceptorClasses.
%
% The last two arguments are the oxygenation fraction and the
% vessel thickness. We set them to be empty here, prompting the
% user to enter these values later.
oxygenationFraction = [];
vesselThickness = [];
fractionBleached = [];
T_receptors = GetHumanPhotoreceptorSS(S, photoreceptorClasses, p.Results.fieldSizeDegrees, p.Results.observerAgeInYears, p.Results.pupilDiameterMm, [], fractionBleached, oxygenationFraction, vesselThickness);

%% Let user choose a photoreceptor class to target

% Here we want to obtain desired contrasts on two
% photoreceptor classes.  We don't have code in place to
% maximize in a specified direction, but we can explicitly
% specify the desired contrasts on the targeted directions.
% One can then maximize by hand, finding out how far one
% can go before it is out of gamut.  A little loop and
% check could automate this process.
whichDirection = 'LMinusM';
whichReceptorsToTarget = [1 2];
whichReceptorsToIgnore = [4 5 6 7 8];
whichReceptorsToMinimize = [];
LMinusMTargetContrast = 0.065;
desiredContrast = [LMinusMTargetContrast -LMinusMTargetContrast];

whichDirection = 'LMS';
whichReceptorsToTarget = [1 2 3];
whichReceptorsToIgnore = [4 5 6 7 8];
whichReceptorsToMinimize = [];
desiredContrast = [1 1 1];

whichDirection = 'S';
whichReceptorsToTarget = [3];
whichReceptorsToIgnore = [4 5 6 7 8];
whichReceptorsToMinimize = [];
desiredContrast = [1];

% User chooses whether to maximize contrast in targeted receptor classes or
% or get it as close to a specified value as possible.
%
% If we target, here we specify the same contrast for all targeted classes.
% This is not necessary, they can differ.  It just makes the demo code a
% bit simpler to yoke them since we only have to prompt for one number.

% Nice message for user
fprintf('\nGenerating stimuli which isolate receptor classes');
for i = 1:length(whichReceptorsToTarget)
    fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToTarget(i)});
end
fprintf('\nGenerating stimuli which ignore receptor classes');
if (~length(whichReceptorsToIgnore) == 0)
    for i = 1:length(whichReceptorsToIgnore)
        fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToIgnore(i)});
    end
else
    fprintf('\n  - None');
end
fprintf('\nThe remaining classes will be silenced\n');


% Don't pin any primaries.  
whichPrimariesToPin = [];

% No smoothness constraint envforced here.  It really wouldn't make
% to much sense for a three-primary monitor, as the smoothness of a
% monitor spectrum is pretty much determined by the spectral shape
% of its primarites.
maxPowerDiff = 10000;



%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    p.Results.primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

%% Compute the contrasts that we got.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);

% Positive modulation of receptors
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Positive modulation contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

% Negative modulation of receptors
modulationReceptors = T_receptors*B_primary*(-(modulationPrimary - backgroundPrimary));
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Negative modulation contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

%% Plots and/or validation check
if ~isdir(p.Results.plotDir)
    mkdir(p.Results.plotDir);
end
curDir = pwd;
cd(p.Results.plotDir);

% Photoreceptor sensitivities
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors,'LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('Sensitivity');
title('Normalized photoreceptor sensitivities');
saveas(theFig1,sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,whichDirection),'pdf');

% Modulation spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);
saveas(theFig2,sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,whichDirection),'pdf');

% Primaries
theFig3 = figure; hold on
plot(modulationPrimary,'r','LineWidth',2);
plot(backgroundPrimary+(-(modulationPrimary-backgroundPrimary)),'g','LineWidth',2);
plot(backgroundPrimary,'k','LineWidth',2);
title('Primary settings');
xlim([0 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary Number (nominal)');
ylabel('Setting');
legend({'Positive', 'Negative', 'Background'},'Location','NorthEastOutside');
saveas(theFig3,sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,whichDirection),'pdf');



%% Return to the directory from whence we started
cd(curDir);
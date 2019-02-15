function resultSet = designNominalSPDs(varargin)
% Nominal primaries and SPDs for isolating post-receptoral mechanisms
%
% Syntax:
%	resultSet = designNominalSPDs
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
%	resultSet             - Cell array of structs. The primaries and SPDs.
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


%% Parse input
p = inputParser;
p.addParameter('calFilePath',getpref('vepMELAanalysis','calFilePath'),@ischar);
p.addParameter('saveDir','~/Desktop/nominalSPDs',@ischar);
p.addParameter('primaryHeadRoom',0,@isscalar)
p.addParameter('observerAgeInYears',25,@isscalar)
p.addParameter('fieldSizeDegrees',30,@isscalar)
p.addParameter('pupilDiameterMm',2,@isscalar)
p.parse(varargin{:});



whichModel = 'human';
whichPrimaries = 'monitor';
curDir = pwd;


% Obtain the cal file and SPDs of the primaries
[calPath, calFileName, ~] = fileparts(p.Results.calFilePath);
cal = LoadCalFile(calFileName,[],calPath);
S = cal.rawData.S;
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

resultSet.photoreceptorClasses = {'L','M','S','Mel','Rod','Lp','Mp','Sp'};

% Make sensitivities.  The wrapper routine is GetHumanPhotoreceptorSS,
% which is in the ContrastSplatter directory.  Each row of the matrix
% T_receptors provides the spectral sensitivity of the photoreceptor class
% in the corresponding entry of the cell array photoreceptorClasses.
%
% The last two arguments are the oxygenation fraction and the vessel
% thickness. We set them to be empty here.
oxygenationFraction = [];
vesselThickness = [];
fractionBleached = [];
T_receptors = GetHumanPhotoreceptorSS(S, photoreceptorClasses, p.Results.fieldSizeDegrees, p.Results.observerAgeInYears, p.Results.pupilDiameterMm, [], fractionBleached, oxygenationFraction, vesselThickness);

% Obtain the isomerization rate for the receptors by the background
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);

% Store the background properties
resultSet.background.primary = backgroundPrimary;
resultSet.background.spd = B_primary*backgroundPrimary;
resultSet.background.wavelengthsNm = SToWls(S);

% set the receptor sets to isolate
whichDirectionSet = {'L','M','S','LMinusM','LMS'};
whichReceptorsToTargetSet = {[1],[2],[3],[1 2],[1 2 3]};
whichReceptorsToIgnoreSet = {[4 5 6 7 8],[4 5 6 7 8],[4 5 6 7 8],[4 5 6 7 8],[4 5 6 7 8]};
whichReceptorsToMinimizeSet = {[],[],[],[],[]};
desiredContrastSet = {[1],[1],[1],[0.065 -0.065], [1 1 1]};

% Loop over the set of directions for which we will generate modulations
for ss = 1:length(whichDirectionSet)
    
    % Extract values from the cell arrays
    whichDirection = whichDirectionSet{ss};
    whichReceptorsToTarget = whichReceptorsToTargetSet{ss};
    whichReceptorsToIgnore = whichReceptorsToIgnoreSet{ss};
    whichReceptorsToMinimize = whichReceptorsToMinimizeSet{ss};
    desiredContrast = desiredContrastSet{ss};
        
    % Don't pin any primaries.
    whichPrimariesToPin = [];
    
    % No smoothness constraint enforced for the monitor primaries
    maxPowerDiff = 10000;
    
    % Obtain the primary settins for the isolating modulation
    modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
        B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
        p.Results.primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);
    
    % Store the modulation primaries
    resultSet.(whichDirection).modulationPrimary = modulationPrimary;
    
    % Calculate and store the positive and negative receptor contrast
    modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
    contrastReceptors = modulationReceptors ./ backgroundReceptors;
    resultSet.(whichDirection).positiveReceptorContrast = contrastReceptors;
    
    modulationReceptors = T_receptors*B_primary*(-(modulationPrimary - backgroundPrimary));
    contrastReceptors = modulationReceptors ./ backgroundReceptors;
    resultSet.(whichDirection).negativeReceptorContrast = contrastReceptors;
    
    % Calculate and store the spectra
    resultSet.(whichDirection).positiveModulationSPD = B_primary*modulationPrimary;
    resultSet.(whichDirection).negativeModulationSPD = B_primary*(backgroundPrimary-(modulationPrimary - backgroundPrimary));
    resultSet.(whichDirection).wavelengthsNm = SToWls(S);
    
    % Create and save results
    if ~isempty(p.Results.saveDir)
        if ~isdir(p.Results.saveDir)
            mkdir(p.Results.saveDir);
        end
        cd(p.Results.saveDir);
       
        % Create a figure with an appropriate title
        fighandle = figure('Name',whichDirection);
               
        % Modulation spectra
        subplot(2,3,[1 4])
        hold on
        plot(resultSet.(whichDirection).wavelengthsNm,resultSet.(whichDirection).positiveModulationSPD,'k','LineWidth',2);
        plot(resultSet.(whichDirection).wavelengthsNm,resultSet.(whichDirection).negativeModulationSPD,'r','LineWidth',2);
        plot(resultSet.background.wavelengthsNm,resultSet.background.spd,'Color',[0.5 0.5 0.5],'LineWidth',2);
        title('Modulation spectra');
        xlim([380 780]);
        ylim([0 0.01]);
        xlabel('Wavelength');
        ylabel('Power');
        legend({'Positive', 'Negative', 'Background'},'Location','NorthEast');
        
        % Primaries
        subplot(2,3,[2 5])
        c = categorical({'R','G','B'});
        hold on
        plot(c,modulationPrimary,'-kx');
        plot(c,backgroundPrimary+(-(modulationPrimary-backgroundPrimary)),'-rx');
        plot(c,backgroundPrimary,'-x','Color',[0.5 0.5 0.5]);
        title('Primary settings');
        ylim([0 1]);
        xlabel('Primary');
        ylabel('Setting');
        
        % Patch
        subplot(2,3,3)
        rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',modulationPrimary)
        axis square
        axis off
        title('+');
        subplot(2,3,6)
        rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',backgroundPrimary+(-(modulationPrimary-backgroundPrimary)))
        axis square
        axis off
        title('-');
        
        % Save the figure
        saveas(fighandle,sprintf('%s_%s_%s_PrimariesAndSPD.pdf',whichModel,whichPrimaries,whichDirection),'pdf');
    end
end


%% Return to the directory from whence we started
cd(curDir);
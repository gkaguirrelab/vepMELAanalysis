function vepMELAanalysisLocalHook

%  vepMELAanalysisLocalHook
%
% Configure things for working on the  vepMELAanalysis project.
%
% For use with the ToolboxToolbox.
%
% If you 'git clone' vepMELAanalysis into your ToolboxToolbox "projectRoot"
% folder, then run in MATLAB
%   tbUseProject('vepMELAanalysis')
% ToolboxToolbox will set up vepMELAanalysis and its dependencies on
% your machine.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/vepMELAanalysisLocalHook.m
%
% Each time you run tbUseProject('vepMELAanalysis'), ToolboxToolbox will
% execute your local copy of this file to do setup for vepMELAanalysis.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


%% Say hello.
projectName = 'vepMELAanalysis';

%% find the correct python version
if ~strcmp(pyversion,'3.7')
    fprintf('Attempting to change python path to anaconda v3.7...\n')
    try
        pyversion('/Applications/anaconda3/bin/python')
        fprintf('Success\n')
    catch
        warnString = 'Could not find Python anaconda v3.7. Visit https://www.anaconda.com/download/';
        warning('localHook:requiresPythonAnaconda',warnString');
    end
else
    fprintf('Python anaconda v3.7 ready\n')
end

%% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify base paths for materials and data
[~, userID] = system('whoami');
userID = strtrim(userID);
switch userID
    case {'melanopsin' 'pupillab'}
        MELA_dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        MELA_analysisBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/'];
    case {'dhb'}
        MELA_dataBasePath = ['/Users1' '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        MELA_analysisBasePath = ['/Users1/' '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/'];
    otherwise
        materialsBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/TOME_materials/hardwareSpecifications/metropsis/PR670 calibration/'];
        MELA_dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        MELA_analysisBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/'];
        MTRP_dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MTRP_data/'];
end

%% Specify where output goes

if ismac
    % Code to run on Mac plaform
    setpref(projectName,'melaDataPath', MELA_dataBasePath);
    setpref(projectName,'melaAnalysisPath', MELA_analysisBasePath);
    setpref(projectName,'calFilePath', fullfile(materialsBasePath,'MetropsisScreen.mat'));
    setpref(projectName,'mtrpDataPath', MTRP_dataBasePath); 
elseif isunix
    % Code to run on Linux plaform
    setpref(projectName,'melaDataPath', MELA_dataBasePath);
    setpref(projectName,'melaAnalysisPath', MELA_analysisBasePath);
elseif ispc
    % Code to run on Windows platform
    warning('No supported for PC')
else
    disp('What are you using?')
end

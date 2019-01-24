function vepMELAnanalysisLocalHook
% Oonfigure things for working on OneLight projects.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute
%   tbUseProject('vepMELAnanalysis')
% to set up for this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
fprintf('Running vepMELAnanalysis local hook\n');
theApproach = 'vepMELAnanalysis';

%% Define protocols for this approach
theProtocols = DefineProtocolNames;
%% Remove old preferences
if (ispref(theApproach))
    rmpref(theApproach);
end
for pp = 1:length(theProtocols)
    if (ispref(theProtocols{pp}))
        rmpref(theProtocols{pp});
    end
end

%% Specify base paths for materials and data
[~, userID] = system('whoami');
userID = strtrim(userID);
switch userID
    case {'melanopsin' 'pupillab'}
        materialsBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        adminBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_admin/'];
        
    case {'dhb'}
        materialsBasePath = ['/Users1'  '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        dataBasePath = ['/Users1' '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
    case {'nicolas'}
        materialsBasePath = '/Volumes/Manta TM HD/Dropbox (Aguirre-Brainard Lab)/MELA_materials';
        dataBasePath = '/Volumes/Manta TM HD/Dropbox (Aguirre-Brainard Lab)/MELA_data';
    otherwise
        materialsBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        adminBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_admin/'];
        
end

%% Set prefs for data
setpref(theApproach,'DataPath',fullfile(dataBasePath));


%% Prefs for individual protocols
for pp = 1:length(theProtocols)
    % Session record base path
    setpref(theProtocols{pp},'SessionRecordsBasePath',fullfile(getpref(theApproach, 'DataPath'),'Experiments',theApproach,theProtocols{pp},'SessionRecords'));
    
    % Data files base path
    setpref(theProtocols{pp},'DataFilesBasePath',fullfile(getpref(theApproach, 'DataPath'),'Experiments',theApproach,theProtocols{pp},'DataFiles'));
end


clear all, close all %#ok<CLALL>
clc
reset(groot)
tic

% Add paths
addpath('Laws\Rate Laws'); % Rate Laws folder
addpath('Laws\Heat Transfer Laws'); % Rate Laws folder
addpath('matDataX');
addpath('export_fig');

%% Initializing values
% Parameters rarely changed by parameter study
DeltaRctH = -158000;           % reaction enthalpy in J/mol, -16500 PB Paper
numSpeciesInlet = 55;        % number of Species being calculated, needs to be larger than numSpeciesGasPhase
numSpeciesGasPhase = 25;      % numbe r of species considered for properties
rhoWall = 7900;                 % densitiy of the Wall, 316L steel (https://www.fabb-it.de/files/datenblaetter/edelstahl.pdf) in kg/m^3

TRctSwitch = true;              % if true: consider energy balance, else calculate at constant reactor temperature TWall
TWallSwitch = true;            % if true: enthalpy balance of Wall is being calculated
constMIXTURESwitch = true;      % if true: Calculate MIXTURE properties at constant pressure and temperature
kModeSwitch = 4;                % 1: k_i is constant value; 2: k_i is being calculated by heat transfer laws, 3: k_i is being calculated by heat transmittance equation
reportStats = 'on';             % Report statistics solvepde

%% Parameters that may change in parameter study or are dependent of changed
% from here until after meshLength parameters can be copied to
% stationary_exec to transfer case properties
elementSize = 0.008;%0.008999982228728;          % size of a mesh triangle
resTol = 5e-8;                % Residual tolerance of PDE solver
TWall = 273.15 + 100;           % Temperature wall in K
TInlet = 273.15 + 255;          % Temperature inlet in K
pInlet = 30;                    % pressure inlet in bar
kRct = 700;                   % Heat Transfer Coefficient inner Wall in W/m^2/K
alphaRct = 700;                  
kCool = 1e6;                   % Heat Transfer Coefficient outer Wall in W/m^2/K
alphaCool = 1e6;
rhoCat = 117.641;
rateFactor = 1;
aRctFactor = 1;
aWallFactor = 1;
DFactor = 1;

WHSV = 6;                     % factor by which mole flows for WHSV=1 are multiplied to achieve required WHSV
ratioInlet = 2.2;               % Ratio of mole flows of H2 to CO (e.g. 2.1 : 1) at inlet
lengthRct = 0.107;            % length of ractor in m
dRct = 15e-4;                    % inner diameter of reactor in m, slit height in kMode=4
wRct = 42.4115e-3;                   % width of slit reactor in m, only used in kMode=4
volRedRct = 0;%43*1e-9;          % Rct1:130 | volumeRct correction in m^3 for internal structures by substraction, dont forget to set to 0 when not used!
% areaCrossRct = 4.4179e-07;     % inner crosssection area in m^2 % areaCrossRct = pi/4*dRct^2 | dRct = sqrt(4*A/pi)
thicknessWall = 10e-4;              % wall thickness in m
meshFactor = 1;
% meshHeight = 0.008999982228728;              % height of of mesh in m
% meshLength = lengthRct;

% copy until here ^

% current date and time in proper format as string
sDate = datestr(datetime('now'),'yyyy-mm-dd');
sDateTime = datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS');

% Check if there's already result folder for the day, otherwise generate it
if ~isfolder(strcat('results/',char(sDate)))
    mkdir('results/',sDate)
end

%% Input variable parameters and intervals
% Input Parameters "parXPar" and Range(Interval) "parXRan" over which parameter is varied
% caseMode decides how cases for parStudy are created
% 1: n-dimensional array with one case for all possible combinations
% 2: 1-dimensional array where each case corresponds to first, second ...
% entry of the parXRan arrays
caseMode = 2;
% IMPORTANT NOTES:
    % - every "parXPar" needs a corresponding "parXRan" and "iX"!
    % - to not consider variable either empty "parXPar" and "parXRan" 
    %   arrays or delete all 3!
    % - for caseMode 2: parXRan vectors need to be of the same size!
par1Par = 'TInlet';
par1Ran = [498.15 513.15 528.15];
i1 = length(par1Ran);
% par2Par = 'TInlet';
% par2Ran = [528.15 528.15 498.15 498.15 528.15 498.15 528.15 498.15 498.15 528.15 498.15 528.15 528.15 498.15 498.15 528.15 513.15 ];
% i2 = length(par2Ran);
% par3Par = 'WHSV';
% par3Ran = [12 12 6 12 6 12 6 6 12 12 6 6 6 12 6 12 9];
% i3 = length(par3Ran);
% par4Par = 'ratioInlet';
% par4Ran = [2.2 1.7 1.7 1.7 2.2 1.7 1.7 2.2 2.2 1.7 1.7 1.7 2.2 2.2 2.2 2.2 1.95];
% i4 = length(par4Ran);
% parXPar = '';
% parXRan = [];
% iX = length(parXRan);

% Read variable names from workspace and create string from all of them
A = whos;
A = [A.name];
% Find at which place a letter sequence starts with "par", followed by any
% single digit, followed by either "Par" or "Ran".
B = regexp(A,'\par\d+[ParRan]+');

parVars = cell(1,length(B));
for k = 1:length(B)
    % Save the 6 symbols following the starting place 
    parVars{k} = A(B(k):B(k)+6);
    % odd elements are replaced by the variable parameters name
    if rem(k,2) ~= 0
        parVars{k} = eval(parVars{k});
    end
    % even elments are replaced by the variable parameters interval
    if rem(k,2) == 0
        parVars{k} = eval(parVars{k});
    end
end
clear A B

% delete empty cells from parVars
parVars(cellfun(@isempty,parVars)) = [];

% Total number of varied parameters
numVar = length(parVars)/2;

%% Create parameter study object
parStudyVarX = parStudy('dateStamp',sDate,...
                        'timeStamp', sDateTime,...
                        'numVar',numVar,...
                        'parVars',parVars,...
                        'caseMode',caseMode);

parStudyVarX.cases;
                    
%% Initialize reactor objects, pass parameters and execute parameter study
% Aus der Parameterstudienklasse "parStudy" wird das Objekt
% "parStudyVarX" erstellt.
% Mit der Methode "setupPS" werden Instanzen von stationary_onDimRct
% in der Eigenschaft "parStudyRct" (cell array) inizialisiert.
% Dann wird mit der Methode "solvePS" mit einer parfor-Schleife
% die Simulation durchgeführt und das Ergebnis in der Eigenschaft 
% "parStudyRct" abgelegt.

parStudyVarX.setupPS('numSpeciesInlet', numSpeciesInlet,...
                        'numSpeciesGasPhase', numSpeciesGasPhase,...
                        'ReportStatistics', reportStats,...
                        'rctEnthalpy', DeltaRctH,...
                        'wallTemperature', TWall,...
                        'elementSize', elementSize,...
                        'lengthRct', lengthRct,...
                        'wRct',wRct,...
                        'dRct',dRct,...
                        'volRedRct',volRedRct,...
                        'rhoCat',rhoCat,...
                        'thicknessWall', thicknessWall,...
                        'meshFactor', meshFactor,...
                        'resTol', resTol,...
                        'TInlet', TInlet,...
                        'pInlet', pInlet,...
                        'WHSV', WHSV,...
                        'ratioInlet', ratioInlet,...
                        'kRct', kRct,...
                        'alphaRct', alphaRct,...
                        'kCool', kCool,...
                        'alphaCool', alphaCool,...
                        'rateFactor', rateFactor,...
                        'aRctFactor', aRctFactor,...
                        'aWallFactor', aWallFactor,...
                        'DFactor', DFactor,...
                        'constMIXTURESwitch', constMIXTURESwitch,...
                        'TRctSwitch', TRctSwitch, ...
                        'kModeSwitch', kModeSwitch,...
                        'TWallSwitch', TWallSwitch);
                        % 'areaCrossRct', areaCrossRct,...
                        % 'meshHeight', meshHeight,...
                        % 'meshLength', meshLength,...
%                         'TConst', TConst,...
%                         'TCool', TCool,...

                        
% start Readout of command window % diaryPath = strcat('results/',parStudyVarX.sDate,'/',parStudyVarX.sDateTime,'_PS','_logCW.txt');
parStudyVarX.diaryPath = strcat('results/',parStudyVarX.sDate,'/',parStudyVarX.sDateTime,'_PS','_logCW.txt');
diary(parStudyVarX.diaryPath);
diary on

parStudyVarX.solvePS;

diary off
logCW = readcell(parStudyVarX.diaryPath,'Delimiter','DontInputAnyRealDelimitersOtherwiseRegexpDoesntWork');
parStudyVarX.caseStat = cell(1,size(parStudyVarX.parVarArrayLine,2));
% read from innermost command (regexp(logCW,'error...) to outside
% index all cells in logCW containing string indicating error, convert with cellfun to logical with 1 for not empty, find with regexp all digits at the end of the strings, convert the digits to array of double, write '0' below each failed case.
parStudyVarX.caseStat(1,cellfun(@str2double,regexp(logCW(~cellfun(@isempty,regexp(logCW,'error while simulation case'))),'\d+$','match'))) = {0};
% fill empty cells with 1 for successfull simulation
parStudyVarX.caseStat(1,cellfun(@isempty,(parStudyVarX.caseStat(1,:)))) = {1};
% transform cell array to logical
parStudyVarX.caseStat = cell2mat(parStudyVarX.caseStat);
% add row for case numbering
parStudyVarX.caseStat(2,:) = 1:size(parStudyVarX.parVarArrayLine,2);

% rctPDE.firstResid = cell2mat(logCM(2,2));

%% evaluate results
evaluatePS(parStudyVarX);

%% Output tables for single value results (e.g. X_CO)
% choose which parameter is constant for the table
% only pairs of parameterand value are permitted (even number of cells)!
% e.g. {'alphaRct','300', 'alphaCool','1e6'}
tableParFix = {'numSpeciesInlet','55'};
% choose which parameters are being written into the table for each case
% any number of stationary_oneDimRct properties is permitted.
tableParVar = {'elementSize','numNodes','resTol','alphaRct','TInlet','WHSV','dRct','thicknessWall','ratioInlet','pInlet','deltaTRctMax','deltaTWallMax','X_H2Mean','X_COMean','X_TotalMean','SCarbon_CO_C1to4','SCarbon_CO_C1to9','SCarbon_CO_Diesel','SCarbon_CO_C21to50','YCarbon_CO_Diesel','STYDiesel','STYHC','MTYDiesel','MTYHC'};
% create output table
parStudyVarX.tableOutput('save',true,'tableParFix',tableParFix,'tableParVar',tableParVar);

%% plot parameter study results
% plotting of parameters over reactor length only works if vectors are of
% the same length.

% Set default figure values
reset(groot);
set(groot,'DefaultFigurePosition',get(0,'screensize').*[1 1 0.5 0.5])
set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesYGrid', 'on', 'DefaultAxesXGrid', 'off')
set(groot,'DefaultAxesGridLineStyle','--')
set(groot,'DefaultAxesFontSize',14)

parStudyVarX.plotInterpPS('save',true,'norm_n',false);
parStudyVarX.showResultFTPS('save',true);
parStudyVarX.showResultSDiesel('save',true);

% save executable file
copyfile('parStudyExecutable.m',strcat('results/',parStudyVarX.sDate,'/',parStudyVarX.sDateTime, '_PS_exec.m'));

% parStudyVarX.computationTime = toc;
% disp(strcat("Solution time is ", num2str(toc), " s."));

% save parameter study
% create folder path and beginning of filename
filename = strcat('results/',parStudyVarX.sDate,'/',parStudyVarX.sDateTime, '_PS');
for l = 2:2:length(parStudyVarX.parVars)
    % add varied parameters to end of filename
    filename = strcat(filename,'_',string(parStudyVarX.parVars(l-1)));
end
filename = strcat(filename,'.mat');
save(filename, 'parStudyVarX');

% Clear workspace
clearvars -except parStudyVarX;

beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
beep
pause(0.2)
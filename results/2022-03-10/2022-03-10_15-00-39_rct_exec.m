clear all, close all %#ok<CLALL>
clc
reset(groot)
%% Test class oneDim_rctWall
% Create instance of class with example parameters
tic

% Add paths
addpath('Laws\Rate Laws'); % Rate Laws folder
addpath('Laws\Heat Transfer Laws'); % Rate Laws folder
addpath('matDataX');

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
resTol = 5e-7;                % Residual tolerance of PDE solver
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
sDateTime = datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS');
sDate = datestr(datetime('now'),'yyyy-mm-dd');

% Check if there's already result folder for the day, otherwise generate it
if ~isfolder(strcat('results/',char(sDate)))
    mkdir('results/',sDate)
end

% Initialize reactor
rctPDE = stationary_oneDimRct('speciesNames', {'H2';'CO';'N2';'CO2';'H2O'; ...
    'CH4';   'C2H6';  'C3H8';  'C4H10'; 'C5H12'; 'C6H14'; 'C7H16'; 'C8H18'; 'C9H20';  'C10H22'; ...
    'C11H24';'C12H26';'C13H28';'C14H30';'C15H32';'C16H34';'C17H36';'C18H38';'C19H40'; 'C20H42'; ...
    'C21H44';'C22H46';'C23H48';'C24H50';'C25H52';'C26H54';'C27H56';'C28H58';'C29H60'; 'C30H62'; ...
    'C31H64';'C32H66';'C33H68';'C34H70';'C35H72';'C36H74';'C37H76';'C38H78';'C39H80'; 'C40H82'; ...
    'C41H84';'C42H86';'C43H88';'C44H90';'C45H92';'C46H94';'C47H96';'C48H98';'C49H100';'C50H102'}, ...
    'numSpeciesInlet', numSpeciesInlet,...
    'numSpeciesGasPhase', numSpeciesGasPhase,...
    'meshFactor', meshFactor,...
    'elementSize',elementSize,...
    'resTol',resTol,...
    'rctEnthalpy', DeltaRctH,...
    'rhoCat', rhoCat,...
    'inletTemperature', TInlet, ...
    'inletPressure', pInlet,...
    'WHSV', WHSV,...
    'ratioInlet',ratioInlet,...
    'wallTemperature', TWall,...
    'kRct', kRct,...
    'kCool', kCool,...
    'alphaRct',alphaRct,...
    'alphaCool',alphaCool,...
    'ReportStatistics', reportStats,...
    'dRct',dRct,...
    'wRct',wRct,...
    'volRedRct',volRedRct,...
    'lengthRct', lengthRct, ...
    'timeStamp', sDateTime,...
    'dateStamp', sDate,...
    'rateFactor', rateFactor,...
    'DFactor', DFactor,...
    'aRctFactor', aRctFactor,...
    'aWallFactor', aWallFactor,...
    'constMIXTURESwitch', constMIXTURESwitch,...
    'TRctSwitch', TRctSwitch, ...
    'kModeSwitch', kModeSwitch,...
    'TWallSwitch', TWallSwitch, ...
    'thicknessWall', thicknessWall );
%     'Velocity', u,...
%     'DiffCoeff', D,...
%     'kinParam', kinParam,...
%     'inletMoleFlow', nInlet, ...
%     'areaCrossRct', areaCrossRct, ...
%     'meshHeight', meshHeight,...
%     'meshLength', meshLength,...
%     'TConst', TConst,...
%     'TCool', TCool,...

% Execute class functions
rctPDE.initModel
rctPDE.domainMesh
rctPDE.setCoeff
rctPDE.setBC
rctPDE.setIC % Only for initial guess - not for transient calculation

% start Readout of command window
diaryPath = strcat('results/',sDate,'/',sDateTime, '_logCW.txt');
diary(diaryPath);
diary on

% Solve model
rctPDE.solve

% logCM = readcell(diaryPath);
% rctPDE.firstResid = cell2mat(logCM(2,2));

% Solve model
rctPDE.evaluate;

% set figure window size
reset(groot);
set(groot,'DefaultFigurePosition',get(0,'screensize').*[1 1 0.5 0.5])
set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesYGrid', 'on', 'DefaultAxesXGrid', 'off')
set(groot,'DefaultAxesGridLineStyle','--')
set(groot,'DefaultAxesFontSize',14)

% Show results via class functions and save plots
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')% Font of mole flow ylabel in plotInterp is not supported and produces a lot of unnecessary warnings
rctPDE.plotInterp('save', true);
rctPDE.showResultFT('save', true);
rctPDE.plotReactor('save', true);
rctPDE.computationTime = toc;
disp(strcat("Solution time is ", num2str(toc), " s."));

% save executable file
copyfile('Executable.m',strcat('results/',rctPDE.sDate,'/',rctPDE.sDateTime, '_rct_exec.m'));

% save case
filename = strcat('results/',rctPDE.sDate,'/',rctPDE.sDateTime, '_rctPDE.mat');
save(filename, 'rctPDE');

diary off

% Clear workspace
clearvars -except rctPDE;

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
% xVector = rctPDE.lengthRct;
% yVector = linspace(0,0.008999982228728,10);
% [gradx,grady] = evaluateGradient(rctPDE.resultPDE.NodalSolution,xVector,yVector);
% 
% gradx = reshape(gradx,size(xVector));
% grady = reshape(grady,size(yVector));
% 
% quiver(xVector,yVector,gradx,grady)
% xlabel('x')
% ylabel('y')

% quick manual plots
% figure(999)
% plot(rctPDE.xq,rctPDE.TRctOut)
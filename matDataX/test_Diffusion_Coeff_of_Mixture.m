clearvars, close all, clc

gasPhase = MIXTURE('names', {'H2','CO','N2','CO2','H2O'});

% Beliebig gewählte zusammensetzung
% diff = gasPhase.diffusivityMix('p',[20 20],'y',[[0.3; 0.3; 0.2; 0.1; 0.1] [0.25; 0.35; 0.2; 0.1; 0.1]],'T',[500 550])
diffBinary = gasPhase.diffusivityBin('p',[20 20],'T',[500 550])

% % Fast nur H2
% diff = gasPhase.diffusivity('y',[0.96; 0.01; 0.01; 0.01; 0.01])
function [Ri, dRHi] = rateLawKwakViaKirsch(R, T, I, y, p_i, conc, kinParam, nSpecies)
%RATELAWTEMPLATE calculates kinetics of example reaction
%   Source: abc
%   Remarks: xyz
%   
%   p_i     partial pressures as line vector of column number 55 in bar
%   Ri      reaction rates as column vector of row number N in unit. e.g. mol/(kg_cat*s)
%   dRHi    reaction enthalpies in J(mol of every species

Ri = zeros(nSpecies, 1);
dRHi = zeros(nSpecies, 1);
dRHi(2) = -160000; % deltaReactionEnthalpy of CO during FTS in J/mol

% Return reaction rates as column vector
Ri = [r_H2;   % reaction rate of H2
      r_CO;   % reaction rate of CO
      0;      % reaction rate of N2
      0;      % reaction rate of CO2
      r_H2O;  % reaction rate of H2O
      r_C];   % 50 reaction rates of 50 alkanes
  
end


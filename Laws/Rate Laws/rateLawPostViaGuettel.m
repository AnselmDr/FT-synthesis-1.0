function r = rateLawPostViaGuettel(T, F, concH2, kinParam)
%RATELAWPOSTVIAGUETTEL is an implementation of a rate law used in a reactor model.
%   The function receives
%       the temperature T in K (double),
%       the pressure p in bar (double),
%       the concentrations of all species in mol/m^3(double vector), %Check
%   	and a set of all kinetic parameters needed for the rate law
%       calculation in kinParam.
%   The function calculates & returns %   the reaction rate of all species
%   in r in mol/m³/s (double vector with same length than y). %Check


%% Calculate the reation rates
% in this section the reaction rates are calculated.

k0 = 3.107*1e10; % m^3 CO/(m^3 Cat*s)
energyA = 120000; % J/mol
k = k0*exp(-energyA/(8.314*T));
r = F*k*concH2;
end


%% References
%   [1] Post, M. F. M.; van't Hoog, A. C.; Minderhoud, J. K.; Sie, S. T.
%   (1989): Diffusion limitations in fischer‐tropsch catalysts. In: AIChE
%   J. 35 (7), S. 1107–1114.
%
%   [2] Guettel, Robert; Turek, Thomas (2009):
%   Comparison of different reactor types for low temperature
%   Fischer–Tropsch synthesis: A simulation study. In: Chemical Engineering
%   Science 64 (5), S. 955–964. DOI: 10.1016/j.ces.2008.10.059.
%
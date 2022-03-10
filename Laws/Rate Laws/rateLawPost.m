function r = rateLawPost(T, concH2, concCO2, kinParam)
%RATELAWPOST is an implementation of a rate law used in a reactor model.
%   The function receives
%       the temperature T in K (double),
%       the pressure p in bar (double),
%       the concentrations of all species in mol/m^3(double vector), %Check
%   	and a set of all kinetic parameters
%   needed for the rate law calculation in kinParam. The function calculates & returns
%   the reaction rate of all species in r in mol/m³/s (double vector with 
%   same length than y).  %Check


    
    %% Calculate the reation rates
    % in this section the reaction rates are calculated.
    
    %r =
    %-kinParam(1)*I*exp(-2.38e4/(8.314*T))*concH2^kinParam(2)*concCO2^kinParam(3);
    %%Change
    r = 0;
end

%%References
%   [1] Post, M. F. M.; van't Hoog, A. C.; Minderhoud, J. K.; Sie, S. T.
%   (1989): Diffusion limitations in fischer‐tropsch catalysts. In: AIChE
%   J. 35 (7), S. 1107–1114.
%
function rates = RateLaw210331(R, T, I, y, p_i, conc, kinParam)
%RATELAWEXAMPLE is an implementation of a rate law used in a reactor model.
%   The function receives the temperature T in K (double), the pressure p in bar
%   (double), the concentrations of all species in mol/m^3(double vector),
%   the adsorption intensity and a set of all kinetic parameters needed for
%   the rate law calculation in kinParam. The function calculates & returns
%   the reaction rate of all species in r in mol/m³/s (double vector with 
%   same length than y). 

    
    %% Calculate the reation rates
    % in this section the reaction rates are calculated.

    
    rRwgs = kinParam(1)*exp(-kinParam(2)/(R*T))*conc(2)*conc(3);
    rForwardMeoh = kinParam(3)*exp(-kinParam(4)/(R*T))*conc(2)*conc(5);
    rReverseMeoh = kinParam(5)*exp(-kinParam(6)/(R*T))*conc(6);
    
    
     % Calculate the reation rates
            rates = ([0;                                      % recation rate of N2 in mol/m³/s
                    -rRwgs-2*rForwardMeoh+2*rReverseMeoh;     % reaction rate of H2 in mol/m³/s
                    -rRwgs;                                   % reaction rate of CO2 in mol/m³/s
                    +rRwgs;                                   % reaction rate of H2O in mol/m³/s
                    +rRwgs-rForwardMeoh+rReverseMeoh;         % reaction rate of CO in mol/m³/s
                    +rForwardMeoh-rReverseMeoh]);             % reaction rate of CH3OH in mol/m³/s
    
end


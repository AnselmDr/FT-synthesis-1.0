function Ri = rateLawKwackViaKirsch(R, T, p_i, nSpecies)
%RATELAWKWACKVIAKIRSCH calculates kinetics of Fischer-Tropsch reaction
%   Source: Kwack, Park,Bae, Ha, Jun (2011): "Development of a kinetic
%   model of the Fischer-Tropsch synthesis reaction with a cobalt-based catalyst"
%   Remarks: 
%   - neglecting olefin readsorption
%   - only reaction constant for chain initiation and propagation temperature-dependent

% Input arguments:
% R    scalar double: universal gas constant
% T    line vector double: T at each location
% p_i  array double: column vector of partial pressurens for each location
%                    in bar

% Output arguments:
% Ri   array double: column vector equals reaction rate for each species;
%                    columns stands for each location in mol/(kg_cat*s)

% function returns nSpecies*numLocations matrix

Ri = zeros(nSpecies, size(p_i,2));

n_c = 50 ; % number of considered alkanes

% Negative p_i produce error by causing complex results (e.g. gamma_h)
% Set negative p_i to 0
p_i(p_i<0) = 0;
    
%**************************************************************************
% parameter input (from regression)
%**************************************************************************
% kinetic parameters at the reference temperature of 518.15 K (245 °C)
% activation energy
 E=[79.9e3      % 1, E_IN, activation energy of chain initiation, [J/mol] 
    99.5e3];    % 2, E_G, activation energy of chain growth, [J/mol]
  
%  k=[1.18e-04                                  % 1, K_H2_ad, [bar-1]
%     0.602918317*exp(-(E(1)/R)*(1/T-1/518.15)) % 2, k_IN, rate constant of chain initiation, Arrhenius equation k=Aexp(-Ea/RT), [mol/kg/s]
%     18.9                                      % 3, k_CH4, rate constant of CH4 formation, [mol/kg/s]
%     2.39                                      % 4, k_Pn, rate constant of formation of paraffin with carbon number n, [mol/kg/s]
%     5.82e-2                                   % 5, k_CO*K_CO_ad, KCOads is the adsorption equilibrium constant of CO and kCO is the dissociation rate constant of CO, [mol/kg/s/bar]
%     3.75e-01*exp(-(E(2)/R)*(1/T-1/518.15))    % 6, kG, rate constant of chain growth, Arrhenius equation k=Aexp(-Ea/RT), [mol/kg/s]
%     1.08];                                    % 7, k_P2, rate constant of ethane formation, [mol/kg/s]
K_H2_ad = 1.18e-04;
k_IN = 0.602918317*exp(-(E(1)/R)*(1./T-1/518.15));
k_CH4 = 18.9;
k_Pn = 3.75;% 3.75 is value of HK, alternative value: 2.39
k_COK_CO_ad = 5.82e-2;
k_G = 3.75e-01*exp(-(E(2)/R)*(1./T-1/518.15));
k_P2 = 1.08;
%**************************************************************************    
gamma_h = (K_H2_ad*p_i(1,:)).^0.5;
beta1 = k_G.*k_IN./k_Pn;
beta2 = k_IN.*gamma_h;
beta3 = k_COK_CO_ad*p_i(2,:);
gamma_CH2 = (-beta2+sqrt(beta2.^2+4*beta1.*beta3))/2./beta1;
gamma1 = k_IN.*gamma_CH2.*gamma_h./(k_G.*gamma_CH2+k_CH4*gamma_h);
A2 = k_G.*gamma_CH2./(k_G.*gamma_CH2 + k_P2*gamma_h);
A_n = k_G.*gamma_CH2./(k_G.*gamma_CH2 + k_Pn*gamma_h);

k = (2:n_c);

% A_sum needs to be calculated for every point in X (value of T)
A_sum = zeros(1,size(p_i,2));
for i = 1:size(p_i,2)
    A_sum(i) = sum(A_n(i).^(k-2));
end

DEN = (1+gamma_h+gamma_CH2+(1+A_sum.*A2).*gamma1).^2;

r_CH4 = k_CH4*gamma1.*gamma_h./DEN;

r_P2 = k_P2*A2.*gamma1.*gamma_h./DEN;

n = (3:n_c)';
r_Pn = k_Pn*A_n.^(n-2).*A2.*gamma1.*gamma_h./DEN;

r_C = [r_CH4; r_P2; r_Pn];
r_Ciso = zeros(n_c-3,1);
r_Co = zeros(n_c-1,1);

i = (1:n_c)';
r_CO = sum(-r_C.*i)+sum(-r_Ciso.*i(4:end))+sum(-r_Co.*i(2:end));
r_H2 = sum(-r_C.*(2*i+1))+sum(-r_Ciso.*(2*i(4:end)+1))+sum(-r_Co.*(2*i(2:end)));
r_H2O = -r_CO;


% Return reaction rates as column vector of row number 55
% 
Ri = [r_H2;   % reaction rate of H2
      r_CO;   % reaction rate of CO
      zeros(1,size(p_i,2));      % reaction rate of N2
      zeros(1,size(p_i,2));      % reaction rate of CO2
      r_H2O;  % reaction rate of H2O
      r_C];   % 50 reaction rates of 50 alkanes

% Without replacing NaN by 0 error "Unsuitable initial guess U0  (default: U0=0)." is produced
Ri(isnan(Ri)) = 0;

% % artificially boost reaction rate of products
% Ri = Ri*1e3;

if ~isreal(Ri)
    disp('Ri is complex!')
end
end


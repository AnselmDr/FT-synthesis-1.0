function Nu_lam = NuPipeLam(Re, Pr, d_h, L)
%NUPIPELAM returns laminar pipe flow nusselt number to function lamNusseltPipe

%Calculate Graetz number
Gz = Re.*Pr.*d_h./L;

Nu1 = 3.66;

Nu2 = 1.615*(Gz).^(1/3);

Nu3 = (2./(1+22*Pr)).^(1/6).*sqrt(Gz);

b = 0.7;

Nu_lam = (Nu1.^3 + b^3 + (Nu2-b).^3 + Nu3.^3).^(1/3);

end


function Nu_turb = NuPipeTurb(Re, Pr, d_h, L)
%NUPIPETURB returns turbulent pipe flow nusselt number to function lamNusseltPipe

chi = 1./(1.8*log(Re)-1.5).^2;

Nu_turb = (chi/8).*Re.*Pr/(1+12.7*sqrt(chi/8).*(Pr.^(2/3)-1)).*(1+(d_h./L).^(2/3));

end


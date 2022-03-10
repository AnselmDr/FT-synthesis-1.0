function Nu = lawNusseltPipe(Re, Pr, d_h, L)
%LAWNUSSELTPIPE returns Nusselt number for pipe flow of arbitrary Re number
%   According to VDI HA G1 (p. 788, equ 29, 30)

% Plausibility checks
if Re < 0
    warning('Re is negative');
end
if Re > 1e9
    warning('Re is > 1e9');
end

gamma = (Re-2300)/(10000-2300);

% Identification flow regime laminar, turbulent or transitional
if Re <= 2300
    Nu = NuPipeLam(Re, Pr, d_h, L);
else
    if Re > 10000
        Nu = NuPipeTurb(Re, Pr, d_h, L);
    else
        Nu = (1-gamma).*NuPipeLam(2300, Pr, d_h, L) + gamma.*NuPipeTurb(10000, Pr, d_h, L);
    end
end

if Nu > 1000
    Nu = 1000;
end

end


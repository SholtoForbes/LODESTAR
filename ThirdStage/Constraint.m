function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0)
mScale = 1;
mfuel_burn = x/mScale;


ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma] = ThirdStageSim(mfuel_burn,k,j,u, phi0, zeta0);
 
Target = [100000-AltF Alt(1)-(min(Alt)+1000) gamma(end)-deg2rad(1)];
% Target = [100000-AltF gamma(end)-deg2rad(1)];


% Target = [100000-AltF (max(q)+.5*q(1))-q(1) gamma(end)-deg2rad(1)];

% Target = (Alt(end) - 100000)^2 -mpayload*10000;
end
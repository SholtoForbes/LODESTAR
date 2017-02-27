function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0)
mScale = 1;
mfuel_burn = x/mScale;


ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma] = ThirdStageSim(mfuel_burn,k,j,u, phi0, zeta0);

Target = [max(q)-60000 160000-AltF];
% Target = (Alt(end) - 100000)^2 -mpayload*10000;
end
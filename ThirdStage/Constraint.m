function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0)
mScale = 1;
mfuel_burn = x/mScale;


ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma] = ThirdStageSim(mfuel_burn,k,j,u, phi0, zeta0);
 
% Target = [100000-AltF Alt(1)-(min(Alt)+4000) gamma(end)-deg2rad(1)];

% if k >= 36000
% Target = [160000-AltF Alt(1)-(min(Alt)+4000) Alt(end)-567000 gamma(end)-deg2rad(1)]; % over 36km this needs to be a bit looser, however it doesnt actually use this. just needs to search it (not sure why).
% else
% Target = [160000-Alt(end) Alt(1)-(min(Alt)+2000) Alt(end)-567000 gamma(end)-deg2rad(1)];
% end

Target = [160000-Alt(end) Alt(1)-(min(Alt)+2000) Alt(end)-567000 gamma(end)-deg2rad(1)];
% Target = [160000-Alt(end) gamma(end)-deg2rad(1)];
% Target = [160000-Alt(end) Alt(end)-567000 gamma(end)-deg2rad(1)];

Target = [160000-Alt(end) Alt(1)-(min(Alt)+1000) Alt(end)-567000 gamma(end)-deg2rad(1)]; % 36km and over

% Target = [100000-AltF Alt(1)-(min(Alt)+10000) gamma(end)-deg2rad(1)];

% Target = [160000-AltF gamma(end)-deg2rad(1)];


% Target = [100000-AltF (max(q)+.5*q(1))-q(1) gamma(end)-deg2rad(1)];

% Target = (Alt(end) - 100000)^2 -mpayload*10000;
end
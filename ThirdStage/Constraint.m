function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0)



ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle] = ThirdStageSim(x,k,j,u, phi0, zeta0);

% ceq = [x(end-1)-Alt(end-1)/1e7];

% Target = [100000-AltF Alt(1)-(min(Alt)+4000) gamma(end)-deg2rad(1)];

% if k >= 36000
% Target = [160000-AltF Alt(1)-(min(Alt)+4000) Alt(end)-567000 gamma(end)-deg2rad(1)]; % over 36km this needs to be a bit looser, however it doesnt actually use this. just needs to search it (not sure why).
% else
% Target = [160000-Alt(end) Alt(1)-(min(Alt)+2000) Alt(end)-567000 gamma(end)-deg2rad(1)];
% end

% Target = [160000-Alt(end) Alt(1)-(min(Alt)+2000) Alt(end)-567000 gamma(end)-deg2rad(1)];
% Target = [160000-Alt(end) gamma(end)-deg2rad(1)];
% Target = [160000-Alt(end) Alt(end)-567000 gamma(end)-deg2rad(1)];

% Alpha
% AoA_max
% Alpha - AoA_max
AoA_constraint = max(Alpha - AoA_max);

Vec_angle_constraint = max(Vec_angle - deg2rad(25))*1e1;

% Target = [100000-Alt(end) max(q)-60000 Alt(end)-567000 gamma(end)-deg2rad(1) AoA_constraint];

Target = [100000-Alt(end) Alt(end)-567000 gamma(end)-deg2rad(1) Vec_angle_constraint];

% Target = [100000-Alt(end) Alt(end)-567000 gamma(end)-deg2rad(1) AoA_constraint];

% Target = [160000-Alt(end) max(q)-60000 Alt(end)-567000 gamma(end)-deg2rad(1)]; 

% Target = [gamma(end)-deg2rad(1) AoA_constraint]; 
%  max(q)-70000    Alt(1)-(min(Alt)) Alt(end)-567000

% ceq = 160000-Alt(end);


% Target = [100000-Alt(end-1) Alt(1)-(min(Alt)) Alt(end-1)-567000 gamma(end-1)-deg2rad(1)];
% x
% Alt
% Target = [100000-Alt(end-1) (max(Alpha)-AoA_max)*10 Alt(1)-(min(Alt)) Alt(end-1)-567000 gamma(end-1)-deg2rad(1)]; % 36km and over
% Target = [100000-Alt(end) -(min(gamma)) Alt(end)-567000 gamma(end)-deg2rad(1)]; % 36km and over

% Target = [160000-Alt(end) (max(q))-80000 Alt(end)-567000 gamma(end)-deg2rad(1)]; 


% Target = [100000-AltF Alt(1)-(min(Alt)+10000) gamma(end)-deg2rad(1)];

% Target = [160000-AltF gamma(end)-deg2rad(1)];


% Target = [100000-AltF (max(q)+.5*q(1))-q(1) gamma(end)-deg2rad(1)];

% Target = (Alt(end) - 100000)^2 -mpayload*10000;
end
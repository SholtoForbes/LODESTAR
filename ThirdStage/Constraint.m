function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0, lb,num_div)



ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div);

% ceq = AltF - 110000;

AoA_constraint = max(Alpha - AoA_max);

Vec_angle_constraint = max(Vec_angle - deg2rad(25))*1e3;


Target = [(100000-Alt(end)) (Alt(end)-400000) gamma(end)-deg2rad(1) Vec_angle_constraint];

end
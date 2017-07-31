function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0, lb,num_div)



ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div);

% ceq = AltF - 110000;

AoA_constraint = max(Alpha - AoA_max);

Vec_angle_constraint = (Vec_angle - deg2rad(20))*1e2;

if length(Vec_angle_constraint) < 150
    Vec_angle_constraint = [Vec_angle_constraint zeros(1,(150-length(Vec_angle_constraint)))]; % in case the rocket crashes and doesnt reach 100s flight
end

Target = [(100000-Alt(end)) (Alt(end)-400000) gamma(end)-deg2rad(1) Vec_angle_constraint(1:150)];
% Target = [(100000-Alt(end)) (Alt(end)-400000) gamma(end)-deg2rad(1)];
end
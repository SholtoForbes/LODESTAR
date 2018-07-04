function [Target,ceq] = Constraint(x,k,j,u, phi0, zeta0, lb,num_div,plotflag)



ceq = 0;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div,plotflag);

% ceq = AltF - 110000;

AoA_constraint = Alpha - AoA_max;
if length(AoA_constraint) < 250
    AoA_constraint = [AoA_constraint zeros(1,(250-length(AoA_constraint)))]; % in case the rocket crashes and doesnt reach 100s flight
end

Vec_angle_constraint = (Vec_angle - deg2rad(8))*1e6;

if length(Vec_angle_constraint) < 250
    Vec_angle_constraint = [Vec_angle_constraint zeros(1,(250-length(Vec_angle_constraint)))]; % in case the rocket crashes and doesnt reach 100s flight
end

AoA_diff_constraint = abs(x(2:end-2)-x(1:end-3))-deg2rad(3); % Constrain the difference beatween two nodes AoA to be less that 3 deg. 

EndofBurnTurn_constraint = x(end)*1000 - x(end-1)*10000/(9.86977*1.5) + 10; % Constrain the AoA variation to end at lease 10s before burn ends. 

Target = [(100000-Alt(end)) (Alt(end)-400000) gamma(end)-deg2rad(1) Vec_angle_constraint(1:250) AoA_constraint(1:250) AoA_diff_constraint EndofBurnTurn_constraint];
% Target = [(100000-Alt(end)) (Alt(end)-400000) gamma(end)-deg2rad(1)];
end
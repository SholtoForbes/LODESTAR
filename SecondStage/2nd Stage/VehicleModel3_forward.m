function ydot = VehicleModel3_forward(f_t, f_y,auxdata,Alpha,Alphadot)

r = f_y(1);
v = f_y(2);
gamma = f_y(3);
m = f_y(4);
phi = f_y(5);
zeta = f_y(6);
% xi = 0;
% phi = 0;
% zeta = deg2rad(97);

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max] = ThirdStageDyn(r,gamma,v,m,Alpha,f_t,auxdata,Alphadot,phi,zeta);

ydot = [rdot;vdot;gammadot;-mdot;phidot;zetadot];
% =========================================================================
end









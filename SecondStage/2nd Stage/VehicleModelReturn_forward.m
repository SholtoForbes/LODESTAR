function ydot = VehicleModelReturn_forward(f_t, f_y,auxdata,alpha,eta,throttle)

% r = f_y(1)+auxdata.Re;
alt = f_y(1);
gamma = f_y(2);
v = f_y(3);
zeta = f_y(4);
phi = f_y(5);
xi = f_y(6);
mFuel = f_y(7);
%

[rdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,0,0,0,1);

ydot = [rdot;gammadot;a;zetadot;phidot;xidot;-Fueldt];
% =========================================================================
end









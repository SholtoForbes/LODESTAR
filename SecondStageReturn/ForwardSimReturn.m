function dz = ForwardSimReturn(y,alphadot,etadot,Atmosphere,interp,flapdeflection,mSPARTAN_empty)

V = y(1);
phi = y(2);
gamma = y(3);
v = y(4);
zeta = y(5);
m = mSPARTAN_empty;
xi = y(6); % longitude doesnt matter
r = V + 6371000;

eta =  y(7);
alpha =  y(8);


A = 62.77; % reference area (m^2)


c = spline( Atmosphere(:,1),  Atmosphere(:,5), V); % Calculate speed of sound using atmospheric data

rho = spline( Atmosphere(:,1),  Atmosphere(:,4), V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

% Calculate Thrust Component ==================================

Thrust = 0;

%======================================================================

Cl1 = interp.Cl_spline(M,alpha);

% body_pitchingmoment = interp.pitchingmoment_spline(M, alpha);% first approximation of pitchingmoment using only body lift

% Flap_lift =q./50000*interp.flaplift_spline(M,alpha,flapdeflection);
Flap_lift = 0;

lift = Cl1*A*q + Flap_lift ;

% Drag = interp.Cd_spline(M,alpha)*A*q +  q/50000*interp.flapdrag_spline(M,alpha,flapdeflection);
Drag = interp.Cd_spline(M,alpha)*A*q;

[rdot,xidot,phidot,gammadot,vdot,zetadot,total_lift] = RotCoordsReturn(r,xi,phi,gamma,v,zeta,lift,Drag,Thrust,m,alpha,eta);

dz = [rdot;phidot;gammadot;vdot;zetadot;xidot;etadot;alphadot];
end
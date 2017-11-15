function [altdot, dfuel, Fueldt, vdot, q, M, D, Thrust, flapdeflection, rho,zeta,phi,eq,zetadot,xi, gammadot] = VehicleModel(time, theta, V, v, mfuel,scattered, const, Atmosphere,zeta,auxdata,alpha,gamma)

% =======================================================
% Vehicle Model for SPARTAN Scramjet Accelerator
% =======================================================
A = auxdata.Stage2.Aref; % reference area (m^2)

eta = deg2rad(0)*ones(1,length(time)); % Roll angle, positive anti-clockwise

dt_array = time(2:end)-time(1:end-1); % Time change between each node pt

m = mfuel + auxdata.Stage2.mStruct + auxdata.Stage3.mTot;


% Aero Data =============================================================
c = ppval(auxdata.interp.c_spline, V); % Calculate speed of sound using atmospheric data

rho = ppval(auxdata.interp.rho_spline, V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

% [Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search);
% [Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search,c,rho,Atmosphere);

Cd_noflaps = auxdata.interp.Cd_spline(M,rad2deg(alpha));
% Cd_noflaps = Cd_noflaps - 0.008; % Compensate for boat tail and engines
Cl_noflaps = auxdata.interp.Cl_spline(M,rad2deg(alpha));
flapdeflection = 0;

D = 0.5*(Cd_noflaps).*A.*rho.*v.^2;
L = 0.5*(Cl_noflaps).*A.*rho.*v.^2;

% Alpha
if const == 14
    D = 1.1*D; % for L/D testing 
end


%Fuel Cost ===========================================================================
% lift = Fl;

T0 = ppval(auxdata.interp.T0_spline, V); 
P0 = ppval(auxdata.interp.P0_spline, V); 

[Isp,Fueldt,eq] = RESTint(M, alpha, auxdata,T0,P0);

for i = 1:length(time)
  if q(i) < 20000
        Isp(i) = Isp(i)*gaussmf(q(i),[1000,20000]);
  end  
end

Thrust = Isp.*Fueldt*9.81; % Thrust in direction of motion

fuelchange_array = -Fueldt(1:end-1).*dt_array ;

dfuel = sum(fuelchange_array); %total change in 'fuel' this is negative



%===================================================
%Rotational Coordinates 
%===================================================

xi = zeros(1,length(time));
phi = zeros(1,length(time));

phi(1) = -0.264;

r = V + 6371000;
i= 1;

% [xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);
[altdot(i),xidot(i),phidot(i),gammadot(i),vdot(i),zetadot(i)] = RotCoordsForward(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),Thrust(i),m(i),alpha(i));
for i = 2:length(time)
xi(i) = xi(i-1) + xidot(i-1)*(time(i) - time(i-1));
phi(i) = phi(i-1) + phidot(i-1)*(time(i) - time(i-1));
% zeta(i) = zeta(i-1) + zetadot(i-1)*(time(i) - time(i-1));

% [xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);

[altdot(i),xidot(i),phidot(i),gammadot(i),vdot(i),zetadot(i)] = RotCoordsForward(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),Thrust(i),m(i),alpha(i));
end

% =========================================================================

end









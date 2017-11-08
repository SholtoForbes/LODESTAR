function [dfuel, Fueldt, a, q, M, Fd, Thrust, flapdeflection, Alpha, rho,lift,zeta,phi,eq,zetadot,xi] = VehicleModel(time, theta, V, v, mfuel,scattered, const,thetadot, Atmosphere,zeta,mstruct,mThirdStage,auxdata)

% =======================================================
% Vehicle Model for SPARTAN Scramjet Accelerator
% =======================================================
A = 62.77; % reference area (m^2)

eta = deg2rad(0)*ones(1,length(time)); % Roll angle, positive anti-clockwise

dt_array = time(2:end)-time(1:end-1); % Time change between each node pt

m = mfuel + mstruct + mThirdStage ;

%===================================================
%Rotational Coordinates 
%===================================================

xi = zeros(1,length(time));
phi = zeros(1,length(time));

phi(1) = -0.264;

r = V + 6371000;
i= 1;

[xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);

for i = 2:length(time)
xi(i) = xi(i-1) + xidot(i-1)*(time(i) - time(i-1));
phi(i) = phi(i-1) + phidot(i-1)*(time(i) - time(i-1));
% zeta(i) = zeta(i-1) + zetadot(i-1)*(time(i) - time(i-1));

[xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);

end

% Aero Data =============================================================
c = ppval(auxdata.interp.c_spline, V); % Calculate speed of sound using atmospheric data

rho = ppval(auxdata.interp.rho_spline, V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

% [Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search);
[Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search,c,rho,Atmosphere);


% Alpha
if const == 14
    Fd = 1.1*Fd; % for L/D testing 
end

% THRUST AND MOTION ==================================================================

%Fuel Cost ===========================================================================
lift = lift_search;

T0 = ppval(auxdata.interp.T0_spline, V); 
P0 = ppval(auxdata.interp.P0_spline, V); 

[Isp,Fueldt,eq] = RESTint(M, Alpha, auxdata,T0,P0);

for i = 1:length(time)
  if q(i) < 20000
        Isp(i) = Isp(i)*gaussmf(q(i),[1000,20000]);
  end  
end

Thrust = Isp.*Fueldt*9.81.*cos(deg2rad(Alpha)); % Thrust in direction of motion

fuelchange_array = -Fueldt(1:end-1).*dt_array ;

dfuel = sum(fuelchange_array); %total change in 'fuel' this is negative

mu_E = 3.986e14; % m^3/s^2 Earth Gravitational Parameter
omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate

a = Thrust./(m) - mu_E.*sin(theta)./r.^2 -Fd./m + omega_E.^2.*r.*cos(phi).*(cos(phi).*cos(theta)+sin(phi).*sin(theta).*sin(zeta));

% =========================================================================

end









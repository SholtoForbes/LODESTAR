function [rdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L] = VehicleModelReturn(gamma, r, v,auxdata,zeta,phi,xi,alpha,eta)

% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.A; % reference area (m^2)

% eta = .0*ones(1,length(time)); % Roll angle

% eta = 0.3 - 0.0001*time;

%Gravity
g = 9.81;

% dt_array = time(2:end)-time(1:end-1); % Time change between each node pt

V = r - auxdata.Re;

m = auxdata.mass;

%===================================================
%
% SECOND STAGE
%
%===================================================

%======================================================
speedOfSound = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,5),V);
mach = v./speedOfSound;
density = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,4),V);
% interpolate coefficients
Cd = auxdata.interp.Cd_spline(mach,rad2deg(alpha));
Cl = auxdata.interp.Cl_spline(mach,rad2deg(alpha));

%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*density.*v.^2;
L = 0.5*Cl.*A.*density.*v.^2;

% D=D*2;
% L=L*3;

%Rotational Coordinates =================================================
%=================================================



i= 1;

T =0;

[rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m,alpha(i),eta(i));

for i = 2:length(r)
[rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m,alpha(i),eta(i));
end

% Aero Data =============================================================
c = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,5), V); % Calculate speed of sound using atmospheric data

rho = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,4), V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

%-heating---------------------------
% kappa = 1.7415e-4;
% Rn = 1; %effective nose radius (m) (need to change this, find actual value)
% 
% heating_rate = kappa*sqrt(rho./Rn).*v.^3; %watts
% 
% Q = zeros(1,length(time));
% Q(1) = 0;
% 
% for i = 1:length(dt_array)
%     Q(i+1) = heating_rate(i)*dt_array(i) + Q(i);
% end



v_H = v.*cos(gamma);

% =========================================================================
end









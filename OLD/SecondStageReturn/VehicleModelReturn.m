function [rdot,xi,phi,gammadot,a,zeta, q, M, D, rho,L] = VehicleModelReturn(time, gamma, V, v, nodes,interp, Atmosphere,zeta,phi,xi,alpha,eta,throttle)

% =======================================================
% Vehicle Model
% =======================================================
A = 62.77; % reference area (m^2)

% eta = .0*ones(1,length(time)); % Roll angle

% eta = 0.3 - 0.0001*time;

%Gravity
g = 9.81;

% dt_array = time(2:end)-time(1:end-1); % Time change between each node pt


mstruct = 4910.5 - 132.8 + 179.41; % mass of everything but fuel from dawids work

m = mstruct;

%===================================================
%
% SECOND STAGE
%
%===================================================

%======================================================
speedOfSound = spline(Atmosphere(:,1),Atmosphere(:,5),V);
mach = v./speedOfSound;
density = spline(Atmosphere(:,1),Atmosphere(:,4),V);
% interpolate coefficients
Cd = interp.Cd_spline(mach,rad2deg(alpha));
Cl = interp.Cl_spline(mach,rad2deg(alpha));

%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*density.*v.^2;
L = 0.5*Cl.*A.*density.*v.^2;

% D=D*2;
% L=L*3;



%Rotational Coordinates =================================================
%=================================================


r = V + 6371000;
i= 1;

T =0;

[rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m,alpha(i),eta(i));

for i = 2:length(time)
xi(i) = xi(i-1) + xidot(i-1)*(time(i) - time(i-1));
phi(i) = phi(i-1) + phidot(i-1)*(time(i) - time(i-1));
zeta(i) = zeta(i-1) + zetadot(i-1)*(time(i) - time(i-1));

[rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m,alpha(i),eta(i));
end

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








function ydot = VehicleModelReturn_forward(f_t, f_y,auxdata,alpha,eta,throttle,mFuelinit,mFuelend)

alt = f_y(1);
gamma = f_y(2);
v = f_y(3);
zeta = f_y(4);
phi = f_y(5);
xi = f_y(6);
mFuel = f_y(7);
% 
% % =======================================================
% % Vehicle Model
% % =======================================================
% A = 62.77; % reference area (m^2)
% 
% % eta = .0*ones(1,length(time)); % Roll angle
% 
% % eta = 0.3 - 0.0001*time;
% 
% %Gravity
% g = 9.81;
% 
% % dt_array = time(2:end)-time(1:end-1); % Time change between each node pt
% 
% 
% mstruct = 4910.5 - 132.8 + 179.41; % mass of everything but fuel from dawids work
% 
% m = mstruct;
% 
% %===================================================
% %
% % SECOND STAGE
% %
% %===================================================
% 
% %======================================================
% speedOfSound = spline(Atmosphere(:,1),Atmosphere(:,5),V);
% mach = v./speedOfSound;
% density = spline(Atmosphere(:,1),Atmosphere(:,4),V);
% % interpolate coefficients
% Cd = interp.Cd_spline(mach,rad2deg(alpha));
% Cl = interp.Cl_spline(mach,rad2deg(alpha));
% 
% %%%% Compute the drag and lift:
% 
% D = 0.5*Cd.*A.*density.*v.^2;
% L = 0.5*Cl.*A.*density.*v.^2;
% 
% % D=D*2;
% % L=L*3;
% 
% %Rotational Coordinates =================================================
% %=================================================
% 
% 
% r = V + 6371000;
% i= 1;
% 
% T =0;
% 
% [rdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(r,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);
% 
% 
% % Aero Data =============================================================
% c = spline( Atmosphere(:,1),  Atmosphere(:,5), V); % Calculate speed of sound using atmospheric data
% 
% rho = spline( Atmosphere(:,1),  Atmosphere(:,4), V); % Calculate density using atmospheric data
% 
% q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure
% 
% M = v./c; % Calculating Mach No (Descaled)
% 
% %-heating---------------------------
% % kappa = 1.7415e-4;
% % Rn = 1; %effective nose radius (m) (need to change this, find actual value)
% % 
% % heating_rate = kappa*sqrt(rho./Rn).*v.^3; %watts
% % 
% % Q = zeros(1,length(time));
% % Q(1) = 0;
% % 
% % for i = 1:length(dt_array)
% %     Q(i+1) = heating_rate(i)*dt_array(i) + Q(i);
% % end
% 
% v_H = v.*cos(gamma);

[rdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,1,1);

ydot = [rdot;gammadot;a;zetadot;phidot;xidot;-Fueldt];
% =========================================================================
end









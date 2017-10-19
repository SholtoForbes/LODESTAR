function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T] = VehicleModelReturn(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,ThirdStage)

interp = auxdata.interp;
% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.A; % reference area (m^2)

% eta = .0*ones(1,length(time)); % Roll angle

% eta = 0.3 - 0.0001*time;

%Gravity
g = 9.81;

% dt_array = time(2:end)-time(1:end-1); % Time change between each node pt

if ThirdStage == 1
m = auxdata.Stage2.mStruct+mFuel+auxdata.Stage3.mTot; 
else
m = auxdata.Stage2.mStruct+mFuel;
end


%===================================================
%
% SECOND STAGE
%
%===================================================


%======================================================

% speedOfSound = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,5),alt);

% density = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,4),alt);






% D=D*2;
% L=L*3;

%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
mach = v./c;
rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

T0 = ppval(interp.T0_spline, alt); 

P0 = ppval(interp.P0_spline, alt);

%% Aerodynamics
% interpolate coefficients
Cd = auxdata.interp.Cd_spline(mach,rad2deg(alpha));
Cl = auxdata.interp.Cl_spline(mach,rad2deg(alpha));

%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*rho.*v.^2;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

[Isp,Fueldt,eq] = RESTM12int(M, alpha, auxdata,T0,P0);

% for i = 1:length(r)
%   if q(i) < 20000
%         Isp(i) = Isp(i)*gaussmf(q(i),[1000,20000]);
%   end  
% end
Isp(q<20000) = Isp(q<20000).*gaussmf(q(q<20000),[1000,20000]);
Fueldt(M<5.1) = 0;

Fueldt = Fueldt.*throttle;

T = Isp.*Fueldt*9.81.*cos(deg2rad(alpha)); % Thrust in direction of motion

% fuelchange_array = -Fueldt(1:end-1).*dt_array ;
% 
% dfuel = sum(fuelchange_array); %total change in 'fuel' this is negative


%Rotational Coordinates =================================================
%=================================================



% i= 1;
% 
% [rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m,alpha(i),eta(i));
% 
% for i = 2:length(r)
% [rdot(i),xidot(i),phidot(i),gammadot(i),a(i),zetadot(i)] = RotCoordsReturn(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m,alpha(i),eta(i));
% end

[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);

% Aero Data =============================================================
% c = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data
% 
% rho = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,4), alt); % Calculate density using atmospheric data

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









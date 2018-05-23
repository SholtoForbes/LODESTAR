function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp,flap_deflection] = VehicleModelAscent(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,ThirdStage)

interp = auxdata.interp;
% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.A; % reference area (m^2)

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

% Cd = auxdata.interp.Cd_spline_EngineOn(mach,rad2deg(alpha));
% Cl = auxdata.interp.Cl_spline_EngineOn(mach,rad2deg(alpha)); 

Cd = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cd_spline_EngineOn.noFuel(mach,rad2deg(alpha));
Cl = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cl_spline_EngineOn.noFuel(mach,rad2deg(alpha));
flap_deflection = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.flap_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.flap_spline_EngineOn.noFuel(mach,rad2deg(alpha));

%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*rho.*v.^2;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

[Isp,Fueldt,eq] = RESTint(M, rad2deg(alpha), auxdata,T0,P0);

Isp(q<20000) = Isp(q<20000).*gaussmf(q(q<20000),[1000,20000]);
Fueldt(M<5.0) = 0;

Fueldt = Fueldt.*throttle;

T = Isp.*Fueldt*9.81.*cos(alpha); % Thrust in direction of motion

%Rotational Coordinates =================================================
%=================================================

[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)



v_H = v.*cos(gamma);

% =========================================================================
end









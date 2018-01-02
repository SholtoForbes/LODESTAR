function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp,q1,flap_deflection] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,ThirdStage)

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

% if ThirdStage == 1
% Cd = auxdata.interp.Cd_spline_EngineOn(mach,rad2deg(alpha));
% Cl = auxdata.interp.Cl_spline_EngineOn(mach,rad2deg(alpha));
% else
% % Cd = auxdata.interp.Cd_spline_EngineOff(mach,rad2deg(alpha));
% % Cl = auxdata.interp.Cl_spline_EngineOff(mach,rad2deg(alpha));   
% 
% Cd = (1-gaussmf(throttle,[.05,1])).*auxdata.interp.Cd_spline_EngineOff(mach,rad2deg(alpha)) + gaussmf(throttle,[.05,1]).*auxdata.interp.Cd_spline_EngineOn(mach,rad2deg(alpha));
% Cl = (1-gaussmf(throttle,[.05,1])).*auxdata.interp.Cl_spline_EngineOff(mach,rad2deg(alpha)) + gaussmf(throttle,[.05,1]).*auxdata.interp.Cl_spline_EngineOn(mach,rad2deg(alpha));   
% 
% % Cd = throttle.*auxdata.interp.Cd_spline_EngineOff(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cd_spline_EngineOn(mach,rad2deg(alpha));
% % Cl = throttle.*auxdata.interp.Cl_spline_EngineOff(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cl_spline_EngineOn(mach,rad2deg(alpha));  
% end

% Cd = auxdata.interp.Cd_spline(mach,rad2deg(alpha));
% Cl = auxdata.interp.Cl_spline(mach,rad2deg(alpha)); 

if ThirdStage == 0
    throttle(M<5.0) =   0; % remove throttle points below operable range on return flight
end

if ThirdStage == 1
    % Interpolate between centre of gravity conditions for FUll and Empty
    % fuel, as fuel depletes
    Cd = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cd_spline_EngineOn.noFuel(mach,rad2deg(alpha));
    Cl = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cl_spline_EngineOn.noFuel(mach,rad2deg(alpha));
    flap_deflection = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.flap_spline_EngineOn.fullFuel(mach,rad2deg(alpha)) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.flap_spline_EngineOn.noFuel(mach,rad2deg(alpha));
else  
    %Interpolate between engine on and engine off cases as throttle is
    %adjusted
    Cd = (1-throttle).*auxdata.interp.Cd_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cd_spline_EngineOn.noThirdStage(mach,rad2deg(alpha));
    Cl = (1-throttle).*auxdata.interp.Cl_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cl_spline_EngineOn.noThirdStage(mach,rad2deg(alpha));  
    flap_deflection = (1-throttle).*auxdata.interp.flap_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.flap_spline_EngineOn.noThirdStage(mach,rad2deg(alpha));  
end


%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*rho.*v.^2;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

[Isp,Fueldt,eq,q1] = RESTint(M, alpha, auxdata,T0,P0);
% 
% Isp(q<20000) = Isp(q<20000).*gaussmf(q(q<20000),[1000,20000]);
% Fueldt(M<5.0) = 0;

Isp(q1<20000) = Isp(q1<20000).*gaussmf(q1(q1<20000),[1000,20000]); % rapidly reduce ISP to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
Fueldt(M<5.0) = 0;
Isp(M<5.0) = 0;

Fueldt = Fueldt.*throttle;
% Fueldt = Fueldt.*gaussmf(throttle,[.05,1]);

T = Isp.*Fueldt*9.81.*cos(deg2rad(alpha)); % Thrust in direction of motion

%Rotational Coordinates =================================================
%=================================================

[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)



v_H = v.*cos(gamma);

% =========================================================================
end









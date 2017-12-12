function [rdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,q1] = VehicleModelReturn(gamma, r, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel)
% function [rdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,trim_constraint] = VehicleModelReturn(gamma, r, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,flapdeflection)
% 
interp = auxdata.interp;
% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.Stage2.Aref; % reference area (m^2)

% eta = .0*ones(1,length(time)); % Roll angle

% eta = 0.3 - 0.0001*time;

%Gravity
g = 9.81;

% dt_array = time(2:end)-time(1:end-1); % Time change between each node pt

alt = r - auxdata.Re;

m = auxdata.Stage2.mStruct+mFuel;



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

% Cd_noflaps = auxdata.interp.Cd_spline_EngineOff(mach,rad2deg(alpha));
% Cl_noflaps = auxdata.interp.Cl_spline_EngineOff(mach,rad2deg(alpha));
% Cm_noflaps = auxdata.interp.Cm_spline_EngineOff(mach,rad2deg(alpha));
  

Cd_noflaps = (1-gaussmf(throttle,[.05,1])).*auxdata.interp.Cd_spline_EngineOff(mach,rad2deg(alpha)) + gaussmf(throttle,[.05,1]).*auxdata.interp.Cd_spline_EngineOn(mach,rad2deg(alpha));
Cl_noflaps = (1-gaussmf(throttle,[.05,1])).*auxdata.interp.Cl_spline_EngineOff(mach,rad2deg(alpha)) + gaussmf(throttle,[.05,1]).*auxdata.interp.Cl_spline_EngineOn(mach,rad2deg(alpha));   


if auxdata.const ==2
    Cd_noflaps = Cd_noflaps*1.1;
end
if auxdata.const ==3
    Cd_noflaps = Cd_noflaps*0.9;
end

if auxdata.const ==6
    Cl_noflaps = Cl_noflaps*1.1;
end
if auxdata.const ==7
    Cl_noflaps = Cl_noflaps*0.9;
end
% Cd_AoA0 = auxdata.interp.Cd_spline(mach,0*ones(length(mach),1));
% Cl_AoA0 = auxdata.interp.Cl_spline(mach,0*ones(length(mach),1));
% Cm_AoA0 = auxdata.interp.Cm_spline(mach,0*ones(length(mach),1));
% % 
% % Cl_AoA0_withflaps = interp.flap_momentCl_scattered(mach,-(Cm_noflaps-Cm_AoA0));
% % Cd_AoA0_withflaps = interp.flap_momentCd_scattered(mach,-(Cm_noflaps-Cm_AoA0)) ;
% 
% Flap_deflection = interp.flap_momentdef_scattered(mach,-(Cm_noflaps-Cm_AoA0)) ;
% 
% if Flap_deflection_linear > 20 
%    Cl_AoA0_withflaps = interp.flap_Cl_scattered(mach,20);
%  Cd_AoA0_withflaps = interp.flap_Cd_scattered(mach,20); 
% elseif Flap_deflection_linear < -20 
%     Cl_AoA0_withflaps = interp.flap_Cl_scattered(mach,-20);
%  Cd_AoA0_withflaps = interp.flap_Cd_scattered(mach,-20); 
% else
% Cl_AoA0_withflaps = polyvaln(auxdata.interp.Cl_poly,[mach,-(Cm_noflaps-Cm_AoA0)]);
% Cd_AoA0_withflaps = polyvaln(auxdata.interp.Cd_poly,[mach,-(Cm_noflaps-Cm_AoA0)]);
% end
    

% Flap_deflection = polyvaln(auxdata.interp.flapdeflection_poly,[mach,-(Cm_noflaps-Cm_AoA0)]);

% Cl_AoA0_withflaps = interp.flap_Cl_scattered(mach,flapdeflection);
% Cd_AoA0_withflaps = interp.flap_Cd_scattered(mach,flapdeflection);
% Cm_AoA0_withflaps = interp.flap_Cm_scattered(mach,flapdeflection);

% flap_Cl = Cl_AoA0_withflaps - Cl_AoA0;
% flap_Cd = Cd_AoA0_withflaps - Cd_AoA0;
% flap_Cm = Cm_AoA0_withflaps - Cm_AoA0;

% trim_constraint = flap_Cm - Cm_noflaps;

%%%% Compute the drag and lift:

% D = 0.5*(Cd_noflaps+flap_Cd).*A.*rho.*v.^2;
% L = 0.5*(Cl_noflaps+flap_Cl).*A.*rho.*v.^2;

D = 0.5*(Cd_noflaps).*A.*rho.*v.^2;
L = 0.5*(Cl_noflaps).*A.*rho.*v.^2;

%% Thrust 

[Isp,Fueldt,eq,q1] = RESTint(M, alpha, auxdata,T0,P0);

if auxdata.const == 4
Isp = Isp*1.1;
end
if auxdata.const == 5
Isp = Isp*0.9;
end
% for i = 1:length(r)
%   if q(i) < 20000
%         Isp(i) = Isp(i)*gaussmf(q(i),[1000,20000]);
%   end  
% end


% Isp(q<20000) = Isp(q<20000).*gaussmf(q(q<20000),[1000,20000]);
% Fueldt(M<5.1) = 0;
Isp(q1<20000) = Isp(q1<20000).*gaussmf(q1(q1<20000),[1000,20000]);
Fueldt(M<5.1) = 0;


Fueldt = Fueldt.*throttle;
% Fueldt = Fueldt.*gaussmf(throttle,[.05,1]);

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

[rdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(r,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);

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









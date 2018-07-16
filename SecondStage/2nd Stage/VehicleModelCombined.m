function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp,q1,flap_deflection,heating_rate] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,ThirdStage,forwardflag)

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

if ThirdStage == 0 && forwardflag ==0;
    throttle(M<5.0) =   0; % remove throttle points below operable range on return flight
end

if ThirdStage == 1
    % Interpolate between centre of gravity conditions for FUll and Empty
    % fuel, as fuel depletes
    
%     Cd = auxdata.interp.Cd_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) ;
%     Cl = auxdata.interp.Cl_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) ;
%     flap_deflection = auxdata.interp.flap_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) ;

mFuel_cyltanks = 710;


% CG_withFuel_noThirdStage = 14.518; % CG from CREO with no third stage but full fuel
% CG_cyltanksEmpty_noThirdStage = 14.297;
% CG_cyltanksEmpty_ThirdStage = (CG_cyltanksEmpty_noThirdStage*(4957.1+710)+16.63*3300)/(4957+710+3300); % CG with third stage, but no fuel
% CG_noThirdStage = (CG_withFuel_noThirdStage*(4957+1562)-12.59*1562)/4957;% CG  at no fuel conditition (it is assumed that the fuel for the return is used so as to not change the CG)
% EngineOn_CG_noFuel = (CG_noThirdStage*4957.1+16.63*3300)/(4957+3300); % CG with third stage, but no fuel
% EngineOn_CG_Fuel = (EngineOn_CG_noFuel*8257.1 + 12.59*1562)/(8257.1+1562); % CG with third stage and full fuel
% 
% % Calculate CG variation before and after cylindrical tank depletion (ie,
% % cylindrical tanks are depleted first)
% CG(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks).*(EngineOn_CG_Fuel-CG_cyltanksEmpty_ThirdStage)+CG_cyltanksEmpty_ThirdStage;
% CG(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel)./(auxdata.Stage2.mFuel-mFuel_cyltanks).*(CG_cyltanksEmpty_ThirdStage-EngineOn_CG_noFuel)+EngineOn_CG_noFuel;




Cd(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks).*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks)).*auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);
Cl(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks).*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks)).*auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);
flap_deflection(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks).*auxdata.interp.flap_spline_EngineOn.fullFuel(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(auxdata.Stage2.mFuel-mFuel_cyltanks)).*auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);


Cd(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel_cyltanks-mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)))./(mFuel_cyltanks).*auxdata.interp.Cd_spline_EngineOn.noFuel(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(mFuel_cyltanks)).*auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);
Cl(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel_cyltanks-mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)))./(mFuel_cyltanks).*auxdata.interp.Cl_spline_EngineOn.noFuel(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(mFuel_cyltanks)).*auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);
flap_deflection(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)) = (mFuel_cyltanks-mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)))./(mFuel_cyltanks).*auxdata.interp.flap_spline_EngineOn.noFuel(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000) + (1-(mFuel(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))-mFuel_cyltanks)./(mFuel_cyltanks)).*auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks)),rad2deg(alpha(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))),alt(mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks))/1000);

Cd = Cd.';
Cl = Cl.';
flap_deflection = flap_deflection.';


%     Cd = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cd_spline_EngineOn.noFuel(mach,rad2deg(alpha),alt/1000);
%     Cl = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.Cl_spline_EngineOn.noFuel(mach,rad2deg(alpha),alt/1000);
%     flap_deflection = (mFuel-mFuelend)./(mFuelinit-mFuelend).*auxdata.interp.flap_spline_EngineOn.fullFuel(mach,rad2deg(alpha),alt/1000) + (1-(mFuel-mFuelend)./(mFuelinit-mFuelend)).*auxdata.interp.flap_spline_EngineOn.noFuel(mach,rad2deg(alpha),alt/1000);
else  
    %Interpolate between engine on and engine off cases as throttle is
    %adjusted
%     Cd = (1-throttle).*auxdata.interp.Cd_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cd_spline_EngineOn.noThirdStage(mach,rad2deg(alpha));
%     Cl = (1-throttle).*auxdata.interp.Cl_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.Cl_spline_EngineOn.noThirdStage(mach,rad2deg(alpha)); 
%     flap_deflection = (1-throttle).*auxdata.interp.flap_spline_EngineOff.noThirdStage(mach,rad2deg(alpha)) + throttle.*auxdata.interp.flap_spline_EngineOn.noThirdStage(mach,rad2deg(alpha));  
    
Cd = (1-throttle).*auxdata.interp.Cd_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.Cd_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);
Cl = (1-throttle).*auxdata.interp.Cl_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.Cl_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);  
    flap_deflection = (1-throttle).*auxdata.interp.flap_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.flap_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);  
end

%%%% Compute the drag and lift:
D = 0.5*Cd.*A.*rho.*v.^2*auxdata.dragmod;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

[Isp,Fueldt,eq,q1] = RESTint(M, rad2deg(alpha), auxdata,T0,P0);
% 
% Isp(q<20000) = Isp(q<20000).*gaussmf(q(q<20000),[1000,20000]);
% Fueldt(M<5.0) = 0;

Isp(q1<20000) = Isp(q1<20000).*gaussmf(q1(q1<20000),[1000,20000]); % rapidly reduce ISP to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
% Fueldt(M<5.0) = 0;

% Isp(M<5.1) = 0; % Reduce ISP to 0 at less than Mach 5.1

% Fueldt = Fueldt.*throttle;
Fueldt = Fueldt.*throttle; %

% T = Isp.*Fueldt*9.81.*cos(alpha); % Thrust in direction of motion
T = Isp.*Fueldt*9.81.*cos(alpha).*gaussmf(throttle,[0.1,1]); % Thrust in direction of motion

%Rotational Coordinates =================================================
%=================================================




[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoords(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta,auxdata.delta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)



v_H = v.*cos(gamma);

%%heating---------------------------
% From Conceptual Shape Optimization of Entry Vehicles, Dirkx & Mooj & NASA
% lecture

%using hot wall correction

kappa = 1.7415e-4; % cited as sutton-graves, from nasa lecture
Rn = 0.005; %effective nose radius (m) 

heating_rate = kappa*sqrt(rho./Rn).*v.^3; %W/m^2


% =========================================================================
end









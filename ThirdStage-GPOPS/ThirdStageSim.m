function [AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff] = ThirdStageSim(alt0,gamma0,v0, phi0, xi0, zeta0, m0, auxdata)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos
% alt0,gamma0,v0, phi0, xi0, zeta0, m0
if isnan(alt0)
    alt0 = 80000;
end


% alt0,gamma0,v0, phi0, xi0, zeta0, m0
% Atmosphere = dlmread('atmosphere.txt');
% Aero = dlmread('AeroCoeffs.txt');
% 
% Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));
% 
% Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));
% 
% CP_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,7));
% 
% CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

Atmosphere = auxdata.Atmosphere;
Drag_interp = auxdata.Drag_interp;

Lift_interp = auxdata.Lift_interp;

iteration = 1;


HelioSync_Altitude = 566.89 + 6371; %Same as Dawids

r_E = 6371000; % earth radius

Orbital_Velocity_f = sqrt(398600/(566.89 + 6371))*10^3; %Calculating the necessary orbital velocity with altitude in km

%Reference area
A = 0.95; % diameter of 1.1m
g = 9.806; %standard gravity

%% Define starting condtions
t(1) = 0.;

dt = 1;

i=1;

r(1) = r_E + alt0;

Alt(1) = alt0;

xi(1) = 0;
    
phi(1) = phi0;

gamma(1) = gamma0;

v(1) = v0;

zeta(1) = zeta0;


% the Isp influences the optimal burn mass
% Isp = 437; % from Tom Furgusens Thesis %RL10

Isp = 317.*0.98; %Kestrel, from Falcon 1 users guide, with efficiency reduction

% Isp = 320
% Isp = 446; %HM7B
% Isp = 340; %Aestus 2

%% Define Vehicle Properties
mHS = 130.9; % Heat Shield Mass

% mEng = 100; %RL10
mEng = 52; %Kestrel
% mEng = 165; %HM7B
% mEng = 138; %Aestus 2 / RS72 from https://web.archive.org/web/20141122143945/http://cs.astrium.eads.net/sp/launcher-propulsion/rocket-engines/aestus-rs72-rocket-engine.html

m = m0;

%% Initiate Simulation
exocond = false;

p_spline = spline( Atmosphere(:,1),  Atmosphere(:,3)); % calculate pressure using atmospheric data

c_spline = spline( Atmosphere(:,1),  Atmosphere(:,5)); % Calculate speed of sound using atmospheric data

rho_spline = spline( Atmosphere(:,1),  Atmosphere(:,4)); % Calculate density using atmospheric data

Alpha = 0;

while (gamma(i) >= 0 && t(i) < 2000 || t(i) < 150) && Alt(end) > 20000 
% iterate until trajectory angle drops to 0, as long as 150s has passed (this allows for the trajectory angle to drop at the beginning of the trajectory)
    
        t(i+1) = t(i) + dt;
       



    p(i) = ppval(p_spline, Alt(i));
    
    if Alt(i) < 85000
        c(i) = ppval(c_spline,  Alt(i)); % Calculate speed of sound using atmospheric data
        rho(i) = ppval(rho_spline, Alt(i)); % Calculate density using atmospheric data
        q(i) = 1/2*rho(i)*v(i)^2;
    else
        c(i) = ppval(c_spline, 85000);
        rho(i) = 0;
        q(i) = 0;
    end
    
    


        T = 0;

 

    M(i) = v(i)/c(i);

    CD(i) = Drag_interp(M(i),rad2deg(Alpha));
    
    CL(i) = Lift_interp(M(i),rad2deg(Alpha));
    
%     CN(i) = CN_interp(M(i),rad2deg(Alpha));
%     
%     cP(i) = CP_interp(M(i),rad2deg(Alpha));

    D(i) = 1/2*rho(i)*(v(i)^2)*A*CD(i);
    L(i) = 1/2*rho(i)*(v(i)^2)*A*CL(i); % Aerodynamic lift
%     N(i) = 1/2*rho(i)*(v(i)^2)*A*CN(i);
    
    
    %%
   
    [rdot(i),xidot(i),phidot(i),gammadot(i),vdot(i),zetadot(i)] = RotCoordsRocket(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m,Alpha,0);
    
%     if i == 1 && gammadot < 0 && x(1) ~= 0
%         fprintf('BAD CONDITIONS!');
%     end
   
    if i == 1
        
    r(i+1) = r(i) + rdot(i)*dt;
    
    Alt(i+1) = r(i+1) - r_E;
    
    xi(i+1) = xi(i) + xidot(i)*dt;
    
    phi(i+1) = phi(i) + phidot(i)*dt;
    
    gamma(i+1) = gamma(i) + gammadot(i)*dt;
    
    v(i+1) = v(i) + vdot(i)*dt;
    
    zeta(i+1) = zeta(i) + zetadot(i)*dt;
    
    else
        % Second order taylor's series with forward difference
    r(i+1) = r(i) + rdot(i)*dt + (rdot(i) - rdot(i-1))/2*dt;
    
    Alt(i+1) = r(i+1) - r_E;
    
    xi(i+1) = xi(i) + xidot(i)*dt + (xidot(i) - xidot(i-1))/2*dt;
    
    phi(i+1) = phi(i) + phidot(i)*dt + (phidot(i) - phidot(i-1))/2*dt;
    
    gamma(i+1) = gamma(i) + gammadot(i)*dt + (gammadot(i) - gammadot(i-1))/2*dt;
    
    v(i+1) = v(i) + vdot(i)*dt + (vdot(i) - vdot(i-1))/2*dt;
    
    zeta(i+1) = zeta(i) + zetadot(i)*dt + (zetadot(i) - zetadot(i-1))/2*dt;
    end

    i = i+1;
end

gamma(end) = 0;

AltF = Alt(end);
AltF_actual = Alt(end);
vF = v(end);

if AltF > 566.89*1000
    AltF = 566.89*1000;
    const = gaussmf(AltF,[10000 566.89*1000]);
elseif AltF < 100000
    AltF = 100000;
    const = gaussmf(AltF,[10000 100000]);
else
    const = 1;
end


if exocond == false
%     fprintf('Did not reach exoatmospheric conditions') % enable this if you are doing testing on the third stage
    m(end) = m(end) - mHS;
end

%Hohmann Transfer, from Dawid (3i)

mu = 398600;
Rearth = 6371; %radius of earth

Omega_E = 7.2921e-5 ; % rotation rate of the Earth rad/s

vexo = sqrt((v(end)*sin(zeta(end)))^2 + (v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))^2); %Change coordinate system when exoatmospheric, add velocity component from rotation of the Earth

inc = acos((v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))/vexo);  % initial orbit inclination angle
inc_diff = acos(-((566.89+6371)/12352)^(7/2))-inc;

% v12 = sqrt(mu / (AltF/10^3 + Rearth))*10^3 - vexo + 2*(v(end))*sin(abs(acos(-((566.89+6371)/12352)^(7/2))-inc)/2); % Final term of this is inclination change cost to get into heliosync orbit
v12 = sqrt(mu / (AltF/10^3 + Rearth))*10^3 - vexo; % without inclination change.It is assumed that the third stage will be near the desired inclination.
v23 = sqrt(mu / (AltF/10^3+ Rearth))*(sqrt(2*HelioSync_Altitude/((AltF/10^3 + Rearth)+HelioSync_Altitude))-1)*10^3 ; 

v34 = sqrt(mu / HelioSync_Altitude)*(1 - sqrt(2*(AltF/10^3 + Rearth)/((AltF/10^3 + Rearth)+HelioSync_Altitude)))*10^3;

dvtot = v12 + v23 + v34;

%as this is happening in a vacuum we can compute whole delta v at once for
%fuel usage, tsiolkovsky rocket equation. 

g_standard = 9.806;

m2 = m(end)/(exp(v12/(Isp*g_standard)));

m3 = m2/(exp(v23/(Isp*g_standard)));

m4 = m3/(exp(v34/(Isp*g_standard)));

% mpayload = m4 - (m(1) - mHS - mEng)*0.09 -mEng; % 9% structural mass used, from falcon 1 guide, second stage masses with no fairing
mpayload = m4 - (m(1) - mHS)*0.09; % 9% structural mass used, from falcon 1 guide, second stage masses with no fairing
% mpayload = m4 - (m(1) - mHS)*0.108695 -mEng; % structural mass used from falcon 1 guide, second stage masses, this assumes fariing not included in dry mass

mpayload = mpayload*const;
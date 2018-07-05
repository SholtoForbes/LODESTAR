function [AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,gamma,zeta,phi, inc] = ThirdStageSimPostAt(k,j,u, phi0, zeta0, m)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos

SCALE_Engine = 1; % changes characteristic length

time1 = cputime;

iteration = 1;

        
Starting_Altitude = k;
Starting_Theta = j;

c = [];
CD = [];
CL = [];
M = [];
CN = [];
CA = [];
vx = [];
vy = [];
rho = [];
t= [];
Theta = [];
Alt = [];
mfuel = [];
Hor = [];
D = [];
L = [];
        


HelioSync_Altitude = 566.89 + 6371; %Same as Dawids

r_E = 6371000; % earth radius

Orbital_Velocity_f = sqrt(398600/(566.89 + 6371))*10^3; %Calculating the necessary orbital velocity with altitude in km

%Reference area
A = 0.866; % diameter of 1.05m

g = 9.81; %standard gravity

% the Isp influences the optimal burn mass
% Isp = 437; % from Tom Furgusens Thesis %RL10
Isp = 317; %Kestrel, from Falcon 1 users guide
% Isp = 446; %HM7B
% Isp = 340; %Aestus 2
%% Define starting condtions
t(1) = 0.;

dt_main = 1; %time step
dt = dt_main;

i=1;

r(1) = r_E + k;

Alt(1) = k;

xi(1) = 0;
    
phi(1) = phi0;

gamma(1) = j;

v(1) = u;

zeta(1) = zeta0;

%% Define Vehicle Properties
mHS = 125.65; % Heat Shield Mass

% mEng = 100; %RL10
mEng = 52; %Kestrel
% mEng = 165; %HM7B
% mEng = 138; %Aestus 2 / RS72 from https://web.archive.org/web/20141122143945/http://cs.astrium.eads.net/sp/launcher-propulsion/rocket-engines/aestus-rs72-rocket-engine.html

% mdot = 14.71; %RL10
mdot = 9.86; %Kestrel
% mdot = 14.8105; %HM7B
% mdot = 16.5; %Aestus 2


%% Initiate Simulation

Alpha = 0;

while (gamma(i) >= 0 && t(i) < 2000 || t(i) < 150) && Alt(end) > 20000 
% iterate until trajectory angle drops to 0, as long as 150s has passed (this allows for the trajectory angle to drop at the beginning of the trajectory)

    
        t(i+1) = t(i) + dt;
        
        


    D = 0;
    L = 0; % Aerodynamic lift
    T=0;
    
    %%
   
    [rdot(i),xidot(i),phidot(i),gammadot(i),vdot(i),zetadot(i)] = RotCoordsRocket(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L,D,T,m,Alpha);
    
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
end

%Hohmann Transfer, from Dawid (3i)

mu = 398600;
Rearth = 6371; %radius of earth

Omega_E = 7.2921e-5 ; % rotation rate of the Earth rad/s

vexo = sqrt((v(end)*sin(zeta(end)))^2 + (v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))^2); %Change coordinate system when exoatmospheric, add velocity component from rotation of the Earth

inc = acos((v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))/vexo);  % initial orbit inclination angle

v12 = sqrt(mu / (AltF/10^3 + Rearth))*10^3 - vexo + 2*(v(end))*sin(abs(acos(-((566.89+6371)/12352)^(7/2))-inc)/2); % Final term of this is inclination change cost to get into heliosync orbit

% v23 = sqrt(mu / (AltF/10^3+ Rearth))*(sqrt(2*HelioSync_Altitude/((AltF/10^3 + Rearth)+HelioSync_Altitude))-1)*10^3 + 2*(v(end)+v12)*sin(abs(acos(-((566.89+6371)/12352)^(7/2))-zeta(end))); % Final term of this is inclination change cost to get into heliosynch orbit

v23 = sqrt(mu / (AltF/10^3+ Rearth))*(sqrt(2*HelioSync_Altitude/((AltF/10^3 + Rearth)+HelioSync_Altitude))-1)*10^3 ; 

v34 = sqrt(mu / HelioSync_Altitude)*(1 - sqrt(2*(AltF/10^3 + Rearth)/((AltF/10^3 + Rearth)+HelioSync_Altitude)))*10^3;

dvtot = v12 + v23 + v34;

%as this is happening in a vacuum we can compute whole delta v at once for
%fuel usage, tsiolkovsky rocket equation. 

g_standard = 9.806;

m2 = m(end)/(exp(v12/(Isp*g_standard)));

m3 = m2/(exp(v23/(Isp*g_standard)));

m4 = m3/(exp(v34/(Isp*g_standard)));

mpayload = m4 - (m(1) - mHS)*0.09 -mEng; % 9% structural mass used, from falcon 1 guide, second stage masses with no fairing


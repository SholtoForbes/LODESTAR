function [AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc] = ThirdStageSim(x,k,j,u, phi0, zeta0)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos

% x(end) = x(end)*10000; % de-scale

SCALE_Engine = 1; % changes characteristic length

time1 = cputime;

Atmosphere = dlmread('atmosphere.txt');
Aero = dlmread('AeroCoeffs.txt');

Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));

Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));

CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

iteration = 1;


rho_init = spline( Atmosphere(:,1),  Atmosphere(:,4), k);
c_init = spline( Atmosphere(:,1),  Atmosphere(:,5), k);

q_init = 0.5*rho_init*u^2;
M_init = u/c_init;

%% Set Max AoA
% determine the maximum allowable normal coefficient with a 10 degree limit
% at 50kPa dynamic pressure, and set the max AoA to match this normal
% coefficient
Alt_50 = spline( Atmosphere(:,4),  Atmosphere(:,1), 50000*2/u(1)^2);
c_50 = spline( Atmosphere(:,1),  Atmosphere(:,5), Alt_50);

M_50 = u(1)/c_50;
% M_50 = 2922.8/c_50;% using a constant velocity

CN_50 = CN_interp(M_50,10);
% CN_50 = CN_interp(M_50,5);
AoA_max(1) = deg2rad(Max_AoA_interp(M_init,CN_50*50000/q_init)); %maximum allowable AoA
%%
AoA_init = x(1); 

[k j u];
        
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
        
mfuel_burn = 2600; % this was chosen to last until after atmospheric exit. This does not make a large amount of difference, and is close to the optimal value.
% mfuel_burn = 2500;
mfuel(1) = mfuel_burn;

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

m(1) = 3300;
% m(1) = 3000;

% mdot = 14.71; %RL10
mdot = 9.86; %Kestrel
% mdot = 14.8105; %HM7B
% mdot = 16.5; %Aestus 2

burntime = mfuel_burn/mdot;

%% Initiate Simulation
exocond = false;
Fuel = true;

j = 1;

p_spline = spline( Atmosphere(:,1),  Atmosphere(:,3)); % calculate pressure using atmospheric data

c_spline = spline( Atmosphere(:,1),  Atmosphere(:,5)); % Calculate speed of sound using atmospheric data

rho_spline = spline( Atmosphere(:,1),  Atmosphere(:,4)); % Calculate density using atmospheric data

while (gamma(i) >= 0 && t(i) < 2000 || t(i) < 150) && Alt(end) > 20000 
% iterate until trajectory angle drops to 0, as long as 150s has passed (this allows for the trajectory angle to drop at the beginning of the trajectory)

    if t(i) > x(end)
        dt = 2; % Increase timestep after end of angle of attack variation. Accuracy doesnt matter as much at that point. 
    end

    mfuel_temp = mfuel(i) - mdot*dt;
    
    if mfuel_temp < mdot*dt && Fuel == true
        mdot = mfuel_temp/dt;
        Fuel = false;
    end
    
        t(i+1) = t(i) + dt;
        
%     if t(i) <= x(end) && Fuel == true
%         Alpha(i) = interp1(0:x(end)/(length(x)-2):x(end),x(1:end-1),t(i),'pchip'); % interpolate between angle of attack points using an interior pchip spline
%     else
%         Alpha(i) = 0;
%     end
    
%     elseif t(i) > x(end)
%         Alpha(i) = 0;
%     end
   

    if t(i) <= 200
        Alpha(i) = interp1(0:200/(length(x)-1):200,x(1:end),t(i),'pchip'); % interpolate between angle of attack points using an interior pchip spline
    elseif t(i) > 200 && t(i) <= burntime
        Alpha(i) = -x(end)*(t(i) - 200)/(burntime - 200) + x(end);
    elseif t(i) > burntime
        Alpha(i) = 0;
    end
 
    p(i) = ppval(p_spline, Alt(i));
    
    if Alt(i) < 85000
        c(i) = ppval(c_spline,  Alt(i)); % Calculate speed of sound using atmospheric data
        rho(i) = ppval(rho_spline, Alt(i)); % Calculate density using atmospheric data
    else
        c(i) = ppval(c_spline, 85000);
        rho(i) = 0;
    
    end
    
    q(i) = 1/2*rho(i)*v(i)^2;
    
%     AoA_max(i) = deg2rad(Max_AoA_interp(M_init,CN_50*50000/q(i))); %maximum allowable AoA
    
    if Fuel == true
%         T = Isp*mdot*9.81 + (1400 - p(i))*1.; % Thrust (N), exit pressure from Rocket Propulsion Analysis program.
        T = Isp*mdot*9.81 - p(i)*A; % Thrust (N)
        mfuel(i+1) = mfuel(i) - mdot*dt;
        
        if q(i) < 10 && exocond == false
            m(i+1) = m(i) - mHS - mdot*dt; %release of heat shield, from dawids glasgow paper
            exocond = true;
        else 
            m(i+1) = m(i) - mdot*dt;
        end
        
    else
        T = 0;

        mfuel(i+1) = mfuel(i);
        
        if q(i) < 10 && exocond == false
            m(i+1) = m(i) - mHS; %release of heat shield, from dawids glasgow paper
            exocond = true;
        else 
            m(i+1) = m(i);
        end
    end
 

    M(i) = v(i)/c(i);
    
    %calculate axial and normal coefficient, and then lift and drag
    %coefficients. from Dawid (3i)


    CD(i) = Drag_interp(M(i),rad2deg(Alpha(i)));
    
    CL(i) = Lift_interp(M(i),rad2deg(Alpha(i)));

%     CA(i) = 0.346 + 0.183 - 0.058*M(i)^2 + 0.00382*M(i)^3;
%     
%     CN(i) = (5.006 - 0.519*M(i) + 0.031*M(i)^2)*rad2deg(Alpha(i));
    
    

    D(i) = 1/2*rho(i)*(v(i)^2)*A*CD(i);
    
    L(i) = 1/2*rho(i)*(v(i)^2)*A*CL(i);
   
    [rdot(i),xidot(i),phidot(i),gammadot(i),vdot(i),zetadot(i)] = RotCoordsRocket(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m(i),Alpha(i));
    
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


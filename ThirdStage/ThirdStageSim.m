function [AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x,k,j,u, phi0, zeta0)
x(1) = x(1)*10000;
x(6) = x(6)*1000;
SCALE_Engine = 1; % changes characteristic length

time1 = cputime;

Atmosphere = dlmread('atmosphere.txt');
Aero = dlmread('AeroCoeffs.txt');
%Simulating the Third Stage Rocket Trajectory
Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));

Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));

CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

iteration = 1;


rho_init = spline( Atmosphere(:,1),  Atmosphere(:,4), k);
c_init = spline( Atmosphere(:,1),  Atmosphere(:,5), k);

q_init = 0.5*rho_init*u^2;
M_init = u/c_init;

% CN_50 = 0.5265; %maximum allowable normal force coefficient, (ten degrees AoA at q=50kPa conditions) This is an assumption, to form a baseline for allowable force
% CN_50 = 0.3; %maximum allowable normal force coefficient

%% this determines the maximum allowable normal coefficient with a 10 degree limit at 50kPa dynamic pressure
Alt_50 = spline( Atmosphere(:,4),  Atmosphere(:,1), 50000*2/u(1)^2);
c_50 = spline( Atmosphere(:,1),  Atmosphere(:,5), Alt_50);
M_50 = u(1)/c_50;
CN_50 = CN_interp(M_50,10);
CN_50*50000/q_init;
M_init;
M_50;
AoA_max = deg2rad(Max_AoA_interp(M_init,CN_50*50000/q_init)); %maximum allowable AoA
if rad2deg(AoA_max) < 5
    AoA_max = deg2rad(5);
end

AoA_init = x(2); 

% AoA_init = AoA_max; 
% if AoA_init > deg2rad(20)
%     AoA_init = deg2rad(20); % keep the angle of attack within the set bounds
% elseif AoA_init < deg2rad(0)
% % if AoA_init < deg2rad(0)
%     AoA_init = deg2rad(0);
% end

% if AoA_init > AoA_max
%     AoA_init = AoA_max;
% end
% 
% 
AoA_mid1 = x(3);
% if AoA_mid1 > AoA_max
%     AoA_mid1 = AoA_max;
% end
% 
AoA_mid2 = x(4);
% if AoA_mid2 > AoA_max
%     AoA_mid2 = AoA_max;
% end

% AoA_mid3 = x(5);
% if AoA_mid3 > AoA_max
%     AoA_mid3 = AoA_max;
% end



% AoA_end = x(6); 

AoA_end = x(5); 

% AoA_end = x(2); 
% AoA_end = deg2rad(10);
% if AoA_end > deg2rad(20)
%     AoA_end = deg2rad(20); % keep the angle of attack within the set bounds
% elseif AoA_end < deg2rad(0)
%     AoA_end = deg2rad(0);
% end

% if AoA_end > AoA_max
%     AoA_end = AoA_max;
% end


% AoA_list = [AoA_init AoA_mid1 AoA_end];

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
        


mfuel_burn = x(1);
mfuel(1) = mfuel_burn;

HelioSync_Altitude = 566.89 + 6371; %Same as Dawids

r_E = 6371000; % earth radius

Orbital_Velocity_f = sqrt(398600/(566.89 + 6371))*10^3; %Calculating the necessary orbital velocity with altitude in km

%Reference area
A = 0.866;

g = 9.81; %standard gravity

% Isp = 437; % from Tom Furgusens Thesis %RL10
Isp = 317; %Kestrel
% Isp = 446; %HM7B
%define starting condtions
t(1) = 0.;

dt_main = 1; %time step
dt = dt_main;

i=1;

%

r(1) = r_E + k;

Alt(1) = k;

xi(1) = 0;
    

phi(1) = phi0;

gamma(1) = j;

v(1) = u;

zeta(1) = zeta0;

mHS = 125; % Heat Shield Mass
% mHS = 302;

% mEng = 100; %RL10
mEng = 52; %Kestrel
% mEng = 165; %HM7B

% m(1) = 2400 + 100 + mHS; %vehicle mass, (to match Dawids glasgow paper)
% m(1) = 2400 + mEng + mHS;
% m(1) = 2900;


% m(1) = 3500; 
m(1) = 3300;
% m(1) = 3200;
% mdot = 14.71; %RL10
mdot = 9.8; %Kestrel
% mdot = 14.8105; %HM7B

burntime = mfuel_burn/mdot;


exocond = false;
Fuel = true;

j = 1;


while gamma(i) >= 0 && Alt(end) < 567*1000 || t(i) < 200;
%     if i == 1
    mfuel_temp = mfuel(i) - mdot*dt;
%     else
%         mfuel_temp = mfuel(i-1) - mdot*2*dt;
%     end
    
    if mfuel_temp < mdot*dt && Fuel == true
        
%         dt = mfuel_temp/mdot;
    mdot = mfuel_temp/dt;
        Fuel = false;
        
    end
    
    dt = dt_main;
        t(i+1) = t(i) + dt;
        
%     if t(i) < burntime
% 
% %         Alpha(i) = atan((tan(AoA_list(end)) - tan(AoA_list(1)))/(burntime) * t(i)  + tan(AoA_list(1)));
%     Alpha(i) = spline([0 burntime/3 2*burntime/3 burntime],[AoA_init,AoA_mid1,AoA_mid2,AoA_end],t(i));
%     elseif t(i) >= burntime
%         Alpha(i) = 0;
%     end

%     if exocond == false
% 
%     Alpha(i) = AoA_init;
%     elseif exocond == true
%         Alpha(i) = 0;
%     end
% Alpha(i) = atan((tan(AoA_list(end)) - tan(AoA_list(1)))/(x(4)) * t(i)  + tan(AoA_list(1)));
if t(i) <= x(6)
Alpha(i) = spline([0 x(6)/3 2*x(6)/3 x(6)],[AoA_init,AoA_mid1,AoA_mid2,AoA_end],t(i));
% Alpha(i) = spline([0 x(7)/4 x(7)/2 3*x(7)/4 x(7)],[AoA_init,AoA_mid1,AoA_mid2,AoA_mid3,AoA_end],t(i));
elseif t(i) > x(6)
    Alpha(i) = 0;
end
 
    p(i) = spline( Atmosphere(:,1),  Atmosphere(:,3), Alt(i));
    
    if Alt(i) < 85000
        c(i) = spline( Atmosphere(:,1),  Atmosphere(:,5), Alt(i)); % Calculate speed of sound using atmospheric data

        rho(i) = spline( Atmosphere(:,1),  Atmosphere(:,4), Alt(i)); % Calculate density using atmospheric data
    else
        c(i) = spline( Atmosphere(:,1),  Atmosphere(:,5), 85000); % if altitude is over 85km, set values of atmospheric data not change. i will need to look at this

        rho(i) = 0;
    end
    
    q(i) = 1/2*rho(i)*v(i)^2;
    
    if Fuel == true
        T = Isp*mdot*9.81 + (1400 - p(i))*1.; % Thrust (N), exit pressure from 
        
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
   
    [rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m(i),Alpha(i));
    
%     if i == 1 && gammadot < 0 && x(1) ~= 0
%         fprintf('BAD CONDITIONS!');
%     end
   
%     if i == 1
    r(i+1) = r(i) + rdot*dt;
    
    Alt(i+1) = r(i+1) - r_E;
    
    xi(i+1) = xi(i) + xidot*dt;
    
    phi(i+1) = phi(i) + phidot*dt;
    
    gamma(i+1) = gamma(i) + gammadot*dt;
    
    v(i+1) = v(i) + vdot*dt;
    
    zeta(i+1) = zeta(i) + zetadot*dt;
%     else
%     r(i+1) = r(i-1) + rdot*2*dt;
%     
%     Alt(i+1) = r(i+1) - r_E;
%     
%     xi(i+1) = xi(i-1) + xidot*2*dt;
%     
%     phi(i+1) = phi(i-1) + phidot*2*dt;
%     
%     gamma(i+1) = gamma(i-1) + gammadot*2*dt;
%     
%     v(i+1) = v(i-1) + vdot*2*dt;
%     
%     zeta(i+1) = zeta(i-1) + zetadot*2*dt;
%     
%     end
    i = i+1;
end

AltF = Alt(end);
AltF_actual = Alt(end);


vF = v(end);
mult1 = 1;
if AltF > 566.89*1000
    mult1 = gaussmf(AltF,[50000 566.89*1000]);
    AltF = 566.89*1000;
end
mult2=1;
% if gamma(end) > 0
%     mult2 = 0;
% end

mult3=1;

% if max(q) > 60000
%     mult3 = gaussmf(max(q),[1000 60000]);
% end

% if max(q) > 50000
%     mult3 = gaussmf(max(q),[1000 50000])
% end

%  if min(gamma) < 0
%     mult3 = gaussmf(min(gamma),[.001 0]);
% end

%  if min(Alt) < Alt(1)
%     mult3 = gaussmf(min(Alt),[100 Alt(1)]);
% end

if exocond == false
%     fprintf('Did not reach exoatmospheric conditions')
    m(end) = m(end) - mHS;
end

%Hohmann Transfer, from Dawid (3i)

mu = 398600;
Rearth = 6371; %radius of earth

Omega_E = 7.2921e-5 ; % rotation rate of the Earth rad/s

vexo = sqrt((v(end)*sin(zeta(end)))^2 + (v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))^2); %Change coordinate system when exoatmospheric, add velocity component from rotation of the Earth

inc = asin(v(end)*sin(pi-zeta(end))/vexo);  % orbit inclination angle

v12 = sqrt(mu / (AltF/10^3 + Rearth))*10^3 - vexo;

v23 = sqrt(mu / (AltF/10^3+ Rearth))*(sqrt(2*HelioSync_Altitude/((AltF/10^3 + Rearth)+HelioSync_Altitude))-1)*10^3 + 2*(v(end)+v12)*sin(abs(acos(-((566.89+6371)/12352)^(7/2))-zeta(end))); % Final term of this is inclination change cost

v34 = sqrt(mu / HelioSync_Altitude)*(1 - sqrt(2*(AltF/10^3 + Rearth)/((AltF/10^3 + Rearth)+HelioSync_Altitude)))*10^3;

dvtot = v12 + v23 + v34;

%as this is happening in a vacuum we can compute while delta v at once for
%fuel usage, tsiolkovsky rocket equation. 


g = 9.81;

m2 = m(end)/(exp(v12/(Isp*g)));

m3 = m2/(exp(v23/(Isp*g)));

m4 = m3/(exp(v34/(Isp*g)));

% mpayload = m4 - 247.4 -mEng; % subtract structural mass
% mpayload = m4 - 189 -mEng; % subtract structural mass
mpayload = m4 - (m(1) - mHS)*0.09 -mEng; % 9% structural mass used, from falcon 1 guide, second stage masses with no fairing
% if exocond == false
%     mpayload = 0;
% end
% Alt(end)
% AltF
mult4=1;
% if AltF < 160000
%     mult4 = gaussmf(AltF,[100000 160000]);
% end
% x(1)
% AltF
% q
% mult

% Alt(end)
mpayload = mult1*mult2*mult3*mult4*mpayload;


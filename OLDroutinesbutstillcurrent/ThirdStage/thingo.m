clear all
k = 35500;
j = 0.05;
u = 2840;


SCALE_Engine = 1.0; % changes characteristic length

time1 = cputime;

Atmosphere = dlmread('atmosphere.txt');
Aero = dlmread('AeroCoeffs.txt');
%Simulating the Third Stage Rocket Trajectory
Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));

Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));

Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

iteration = 1;


rho_init = spline( Atmosphere(:,1),  Atmosphere(:,4), k);
c_init = spline( Atmosphere(:,1),  Atmosphere(:,5), k);
q_init = 0.5*rho_init*u^2;
M_init = u/c_init;

CN_50 = 0.5265; %maximum allowable normal force coefficient, (ten degrees AoA at q=50kPa conditions) This is an assumption, to form a baseline for allowable force

AoA_max = deg2rad(Max_AoA_interp(M_init,CN_50))*50000/q_init; %maximum allowable AoA



AoA_list = [deg2rad(16) deg2rad(14)];

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
        


mfuel_burn = 800;
mfuel(1) = mfuel_burn;

HelioSync_Altitude = 566.89 + 6371; %Same as Dawids

r_E = 6371000; % earth radius

Orbital_Velocity_f = sqrt(398600/(566.89 + 6371))*10^3; %Calculating the necessary orbital velocity with altitude in km

% maxturnratepos = deg2rad(8)/80; %rad/s these are from dawids glasgow paper, need to figure out where these come from
% maxturnrateneg = -deg2rad(8)/110; %rad/s
% maxturnrateneg = -0.003;

%Reference area from Dawids Cadin
% A = 0.87;
A = 1.84*SCALE_Engine^2;

g = 9.81; %standard gravity

Isp = 451;

%define starting condtions
t(1) = 0.;

dt_main = .5; %time step
dt = dt_main;

i=1;

%

r(1) = r_E + k;

Alt(1) = k;

xi(1) = 0;
    
phi(1) = -0.1314;

gamma(1) = deg2rad(j);

v(1) = u;

zeta(1) = deg2rad(84);

m(1) = 2700 + 167*SCALE_Engine^3; %vehicle mass, (to match Dawids glasgow paper)

mdot = 22.4*SCALE_Engine^2;

burntime = mfuel_burn/mdot;


exocond = false;
Fuel = true;

j = 1;


while gamma(i) >= 0;
    mfuel_temp = mfuel(i) - mdot*dt;
    if mfuel_temp < mdot*dt && Fuel == true
        
        dt = mfuel_temp/mdot;
        t(i+1) = t(i) + dt;
        Fuel = false;
    else
        dt = dt_main;
        t(i+1) = t(i) + dt;
    end
    
    if t(i) < burntime

        Alpha(i) = atan((tan(AoA_list(end)) - tan(AoA_list(1)))/(burntime) * t(i)  + tan(AoA_list(1)));
    elseif t(i) >= burntime
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
        T = 99100*SCALE_Engine^2;
        
        mfuel(i+1) = mfuel(i) - mdot*dt;
        
        if q(i) < 10 && exocond == false
            m(i+1) = m(i) - 302.8 - mdot*dt; %release of heat shield, from dawids glasgow paper
            exocond = true;
        else 
            m(i+1) = m(i) - mdot*dt;
        end
        
    else
        T = 0;

        mfuel(i+1) = mfuel(i);
        
        if q(i) < 10 && exocond == false
            m(i+1) = m(i) - 302.8; %release of heat shield, from dawids glasgow paper
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

    D(i) = 1/2*rho(i)*(v(i)^2)*A*CD(i);
    
    L(i) = 1/2*rho(i)*(v(i)^2)*A*CL(i);
   
    [rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(r(i),xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T,m(i),Alpha(i));
    
    if i == 1 && gammadot < 0
        fprintf('BAD CONDITIONS!');
    end
   
    
    r(i+1) = r(i) + rdot*dt;
    
    Alt(i+1) = r(i+1) - r_E;
    
    xi(i+1) = xi(i) + xidot*dt;
    
    phi(i+1) = phi(i) + phidot*dt;
    
    gamma(i+1) = gamma(i) + gammadot*dt;
    
    v(i+1) = v(i) + vdot*dt;
    
    zeta(i+1) = zeta(i) + zetadot*dt;
    
    
    i = i+1;
end

AltF = Alt(end);
AltF_actual = Alt(end);


vF = v(end);

if AltF > HelioSync_Altitude
    AltF = HelioSync_Altitude;
end


if exocond == false
%     fprintf('Did not reach exoatmospheric conditions')
    m(end) = m(end) - 302.8;
end

%Hohmann Transfer, from Dawid (3i)

mu = 398600;
Rearth = 6371; %radius of earth

Omega_E = 7.2921e-5 ; % rotation rate of the Earth rad/s

vexo = sqrt((v(end)*sin(zeta(end)))^2 + (v(end)*cos(zeta(end)) + r(end)*Omega_E*cos(phi(end)))^2); %Change coordinate system when exoatmospheric, add velocity component from rotation of the Earth

inc = asin(v(end)*sin(pi-zeta(end))/vexo);  % orbit inclination angle

v12 = sqrt(mu / (AltF/10^3 + Rearth))*10^3 - vexo;

v23 = sqrt(mu / (AltF/10^3+ Rearth))*(sqrt(2*HelioSync_Altitude/((AltF/10^3 + Rearth)+HelioSync_Altitude))-1)*10^3;

v34 = sqrt(mu / HelioSync_Altitude)*(1 - sqrt(2*(AltF/10^3 + Rearth)/((AltF/10^3 + Rearth)+HelioSync_Altitude)))*10^3;

dvtot = v12 + v23 + v34;

%as this is happening in a vacuum we can compute while delta v at once for
%fuel usage, tsiolkovsky rocket equation. this should have gravity maybe


g = 9.81;

m2 = m(end)/(exp(v12/(Isp*g)));

m3 = m2/(exp(v23/(Isp*g)));

m4 = m3/(exp(v34/(Isp*g)));

mpayload = m4 - 247.4 -167*SCALE_Engine^3; % subtract structural mass, from Dawids glasgow paper

if exocond == false
    mpayload = 0;
end

mpayload


function dz = ForwardSim(y,alpha,communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered,liftact,dragact,thrustact)

V = y(1);
phi = y(2);
gamma = y(3);
v = y(4);
zeta = y(5);
m = y(6);
xi = 0; % longitude doesnt matter
r = V + 6371000;

flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,6));


[MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';

Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');

A = 62.77*SPARTAN_SCALE^(2/3); % reference area (m^2)


c = spline( Atmosphere(:,1),  Atmosphere(:,5), V); % Calculate speed of sound using atmospheric data

rho = spline( Atmosphere(:,1),  Atmosphere(:,4), V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

% Calculate Thrust Component ==================================

kpa50_V = interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000./v.^2);% Vitude at each velocity, if it were a constant 50kPa q trajectory (used to compare with communicator matrix results) 
kpa50_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), kpa50_V); % Calculate density using atmospheric data

if const == 1 || const == 14

    temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), V);
    Efficiency = zeros(1,length(q));
    for i = 1:length(q)
        if q(i) < 50000
        Efficiency(i) = rho/(50000*2/v^2); % dont change this
        t_ratio(i) = temp_actual(i)/kpa50_temp(i);
        else
        Efficiency(i) = 1; % for 50kPa
        t_ratio(i) = 1; % 
        end
    end
elseif const == 12
    constq_V = interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000./v.^2); % Vitude at each velocity, if it were a constant q trajectory (used to compare with communicator matrix results) 
    constq_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), constq_V); % Calculate density using atmospheric data
    temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), V);
    Efficiency = zeros(1,length(q));
    for i = 1:length(q)
        if q(i) < 55000
        Efficiency(i) = rho/(50000*2/v^2); % dont change this
        t_ratio(i) = 1;
        else
        Efficiency(i) = 1.1; % for 55kPa
        t_ratio(i) = 1;
        end
    end
elseif const == 13
    constq_V = interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000./v.^2); % Vitude at each velocity, if it were a constant q trajectory (used to compare with communicator matrix results) 
    constq_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), constq_V); % Calculate density using atmospheric data
    temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), V);
    Efficiency = zeros(1,length(q));
    for i = 1:length(q)
        if q(i) < 45000
            Efficiency(i) = rho/(50000*2/v^2); % dont change this
            t_ratio(i) = 1;
        else
            Efficiency(i) = .9; % for 45kPa
            t_ratio(i) = 1;
        end
    end
elseif const == 3 || const == 31
temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), V);
Efficiency = rho./(50000*2./v.^2); % linear rho efficiency, scaled to rho at 50000kpa
t_ratio = 1;
end
Efficiency = 1;
T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), V); 
P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), V); 

[Isp,Fueldt] = RESTM12int(M, alpha, Efficiency, scattered, SPARTAN_SCALE,T0,P0);

Thrust = Isp.*Fueldt*9.81;

%======================================================================

Cl1 = Cl_spline(M,alpha);

body_pitchingmoment = pitchingmoment_spline(M, alpha);% first approximation of pitchingmoment using only body lift

Flap_lift = q./50000*flaplift_spline(M,alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments

% total_lift = Cl1*A*q + Flap_lift + Thrust*sin(deg2rad(alpha)) %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

lift = Cl1*A*q + Flap_lift ;

Drag = Cd_spline(M,alpha)*A*q +  q/50000*flapdrag_spline(M,alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);

flapdeflection = flapdeflection_spline(M,alpha,-body_pitchingmoment);

[rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsForward(r,xi,phi,gamma,v,zeta,lift,Drag,Thrust,m,alpha);

dz = [rdot;phidot;gammadot;vdot;zetadot;-Fueldt];
end
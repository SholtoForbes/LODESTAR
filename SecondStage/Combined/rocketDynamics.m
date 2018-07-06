function [dz,q,xi] = rocketDynamics(z,u,t,phase,interp,Throttle,Vehicle,Atmosphere)
% function [dz,q,phi] = rocketDynamics(z,u,t,phase,scattered)
%This function determines the dynamics of the system as time derivatives,
%from input primals. Called by FirstStageDynamics

global mach

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
alpha = z(5,:);
zeta = z(6,:);

dalphadt = z(7,:);

phi = z(8,:);

dalphadt2 = u(1,:); % control is second derivative of AoA with time

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

if h>=0
density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);
else
    density = interp1(Atmosphere(:,1),Atmosphere(:,4),0);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),0);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),0);
end

q = 0.5*density.*v.^2;

% Merlin 1C engine 
T = Vehicle.T.SL + (101325 - P_atm)*Vehicle.T.Mod; % Thrust from Falcon 1 users guide. exit area calculated in SCALING.docx

T = T.*Throttle; % Throttle down



Isp = Vehicle.Isp.SL + (101325 - P_atm)*Vehicle.Isp.Mod; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide

dm = -T./Isp./g;

mach = v./speedOfSound;

% interpolate coefficients
Cd = interp.DragGridded(mach,rad2deg(alpha));
Cl = interp.LiftGridded(mach,rad2deg(alpha));

%%%% Compute the drag and lift:

Area = Vehicle.Area; 

D = 0.5*Cd.*Area.*density.*v.^2;
global L
L = 0.5*Cl.*Area.*density.*v.^2;

switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    zeta = 0;
    case 'postpitch'
    %Do nothing
end

xi = 0; 

% This determines the dynamics of the system.
% Set up like this because xi is a quasi-forward sim instead of a primal




[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);



switch phase
    case 'prepitch'
    dgamma = 0; % Control Trajectory Angle 
    dzeta = 0;
    case 'postpitch'
    %Do nothing
end

dz = [dr;dv;dm;dgamma;dalphadt;dzeta;dalphadt2;dphi;dxi];

if any(isnan(dz))
    disp('NaN Values Detected')
end

end
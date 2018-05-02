function [dz,q,phi,xi] = rocketDynamics(z,u,t,phase,scattered)

%This function determines the dynamics of the system as time derivatives,
%from input primals. Called by LanderDynamics

global mach
Atmosphere = dlmread('atmosphere.txt');
h = z(1,:);   %Height
v = z(2,:);   %Velocity
gamma = z(3,:);
alpha = z(4,:);
zeta = z(5,:);
dalphadt = z(6,:);

m = 836.58*ones(1,length(h)); % mass (kg)

dalphadt2 = u(1,:); % control is second derivative of AoA with time

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);

q = 0.5*density.*v.^2;

SCALE = 1.;
% SCALE = 1; %this is engine exit area scale
% 5 boosters, engines from Microcosm
% T = 5*(88960 - 0.1404*P_atm);% JC changed
T = 0*ones(1,length(h));
% Patm -> T
% 101325 -> 74730
% 0 -> 88960
% Isp = 285 - 0.0005*P_atm; % JC changed 
% Patm -> Isp
% 101325 -> 239
% 0 -> 285


mach = v./speedOfSound;

% interpolate coefficients
% Cd = scattered.DragGridded(mach,rad2deg(alpha));
% Cl = scattered.LiftGridded(mach,rad2deg(alpha));
Cd = scattered.Drag(mach,rad2deg(alpha));
Cl = scattered.Lift(mach,rad2deg(alpha));

%%%% Compute the drag and lift:
global SPARTANscale
Area = 8.308; %62.77*SPARTANscale^(2/3);  JC changed
D = 0.5*Cd.*Area.*density.*v.^2;
global L
L = 0.5*Cl.*Area.*density.*v.^2;



%%%% Complete the calculation:


% xi = 0*ones(1,length(h)); 
% phi = -0.2138*ones(1,length(h));
% zeta = deg2rad(97)*ones(1,length(h));


switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    zeta = 0;
    case 'postpitch'
    %Do nothing
end

xi = 0; 
phi = -0.2220; % initial latitude

% This determines the dynamics of the system.
% Set up like this because phi is a quasi-forward sim instead of a primal
i = 1;
[dr(i),dxi(i),dphi(i),dgamma(i),dv(i),dzeta(i)] = RotCoords(h(i)+rEarth,xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m(i),alpha(i),phase);
for i= 2:length(t)

xi(i) = xi(i-1) + dxi(i-1)*(t(i) - t(i-1));
phi(i) = phi(i-1) + dphi(i-1)*(t(i) - t(i-1));

[dr(i),dxi(i),dphi(i),dgamma(i),dv(i),dzeta(i)] = RotCoords(h(i)+rEarth,xi,phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m(i),alpha(i),phase);
end


switch phase
    case 'prepitch'
    dgamma = 0; % Control Trajectory Angle 
    dzeta = 0;
    case 'postpitch'
    %Do nothing
end


dz = [dr;dv;dgamma;dalphadt;dzeta;dalphadt2];

if any(isnan(dz))
    disp('NaN Values Detected')
end
% dz = [dr;dv;dm;dgamma];
end
function [dz,q,xi] = rocketDynamics(z,u,t,phase,scattered,Throttle)
% function [dz,q,phi] = rocketDynamics(z,u,t,phase,scattered)
%This function determines the dynamics of the system as time derivatives,
%from input primals. Called by LanderDynamics

global mach
Atmosphere = dlmread('atmosphere.txt');
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

SCALE = 1.;
% SCALE = 1; %this is engine exit area scale
% Merlin 1C engine 
T = 555900*SCALE + (101325 - P_atm)*0.5518*SCALE; % Thrust from Flacon 1 users guide. exit area calculated in SCALING.docx

T = T.*Throttle; % Throttle down

Isp = 275 + (101325 - P_atm)*2.9410e-04; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide

dm = -T./Isp./g*SCALE;


mach = v./speedOfSound;

% interpolate coefficients
Cd = scattered.DragGridded(mach,rad2deg(alpha));
Cl = scattered.LiftGridded(mach,rad2deg(alpha));

%%%% Compute the drag and lift:
global SPARTANscale
Area = 62.77*SPARTANscale^(2/3); 

D = 0.5*Cd.*Area.*density.*v.^2;
global L
L = 0.5*Cl.*Area.*density.*v.^2;

% D=2*D;
% L=2*L;

% ref_length = 30.96;
% Cm = scattered.MomentGridded(mach,rad2deg(alpha));
% moment = 0.5*Cm.*ref_length.*Area.*density.*v.^2;
% thrust_vector = asin(-moment./(8.*T));
% T = T - T.*abs(sin(thrust_vector));
% L = L - T.*sin(thrust_vector); % negative so that nose down pith results in positive lift



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
% phi = -0.2220; % initial latitude

% This determines the dynamics of the system.
% Set up like this because phi is a quasi-forward sim instead of a primal
i = 1;
[dr(i),dxi(i),dphi(i),dgamma(i),dv(i),dzeta(i)] = RotCoords(h(i)+rEarth,xi(i),phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m(i),alpha(i),phase);

for i= 2:length(t)

% phi(i) = phi(i-1) + dphi(i-1)*(t(i) - t(i-1));
xi(i) = xi(i-1) + dxi(i-1)*(t(i) - t(i-1));
[dr(i),dxi(i),dphi(i),dgamma(i),dv(i),dzeta(i)] = RotCoords(h(i)+rEarth,xi,phi(i),gamma(i),v(i),zeta(i),L(i),D(i),T(i),m(i),alpha(i),phase);
end


switch phase
    case 'prepitch'
    dgamma = 0; % Control Trajectory Angle 
    dzeta = 0;
    case 'postpitch'
    %Do nothing
end


dz = [dr;dv;dm;dgamma;dalphadt;dzeta;dalphadt2;dphi];

if any(isnan(dz))
    disp('NaN Values Detected')
end
% dz = [dr;dv;dm;dgamma];
end
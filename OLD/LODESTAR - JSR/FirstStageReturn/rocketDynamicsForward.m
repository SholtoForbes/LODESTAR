function [dz] = rocketDynamicsForward(z,alpha,phase,scattered)
global mach
Atmosphere = dlmread('atmosphere.txt');
h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
zeta = z(5,:);
phi = z(6,:);


dalphadt = 0;

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
% T = 5*(88960 - 0.1404*P_atm);% JC changed
% Patm -> T
% 101325 -> 74730
% 0 -> 88960
% Isp = 285 - 0.0005*P_atm; % JC changed 
% Patm -> Isp
% 101325 -> 239
% 0 -> 285

T = 0*ones(1,length(h));
dm = 0*ones(1,length(h));


mach = v./speedOfSound;
Cd = scattered.Drag(mach,rad2deg(alpha));
Cl = scattered.Lift(mach,rad2deg(alpha));
% Cd = scattered.DragGridded(mach,rad2deg(alpha));
% Cl = scattered.LiftGridded(mach,rad2deg(alpha));

%%%% Compute the drag:
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
    case 'postpitch'
    %Do nothing
end

xi = 0; 


[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoords(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);

% dzeta
if isnan(dgamma)
dgamma = 0;
end



dz = [dr;dv;dm;dgamma;dphi;dzeta];

if any(isnan(dz))
    disp('NaN Values Detected')
end
% dz = [dr;dv;dm;dgamma];
end
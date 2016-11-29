function dz = rocketDynamicsFullSize(z,u,phase,scattered)
global mach

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
alpha = z(5,:);

% alpha = u(1,:);

dalphadt = u(1,:);

if isnan(gamma)
    gamma = 1.5708;
end

% T = u(1,:);        %Thrust
T = 460000;

density = 1.474085291*(0.9998541833.^h);  %Data fit off of wolfram alpha



speedOfSound = 280;  %(m/s)  %At 10 km altitude
mach = v/speedOfSound;

Cd = scattered.Drag(mach,rad2deg(alpha));
Cl = scattered.Lift(mach,rad2deg(alpha));

%%%% Compute the drag:
Area = 62.77;  
D = 0.5*Cd.*Area.*density.*v.^2;
L = 0.5*Cl.*Area.*density.*v.^2;

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

%%%% Complete the calculation:
global Tmax
% dm = -60*ones(1,length(h)).*T/Tmax.*Tmax/200000;   %mass rate
dm = -160*ones(1,length(h)).*T/Tmax;

% alpha = 0*ones(1,length(h));


xi = 0*ones(1,length(h));
phi = 0*ones(1,length(h));
zeta = 0*ones(1,length(h));


switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    case 'postpitch'
    %Do nothing
end


[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoords(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);

if isnan(dgamma)
dgamma = 0;
end

dz = [dr;dv;dm;dgamma;dalphadt];

end
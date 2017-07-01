function XDOT = SecondStageReturnDynamics(primal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global a
global Fueldt
global const
global scale
global zetadot
global nodes
global scattered
global Atmosphere


V = primal.states(1,:)*scale.V ; 

v = primal.states(2,:)*scale.v ; 

vdot = a;
 
gamma = primal.states(3,:)*scale.gamma ; 

% alpha = primal.states(4,:);
alpha = deg2rad(7)*ones(1,length(V));

zeta = primal.states(4,:);

phi = primal.states(5,:);

xi = primal.states(6,:);

eta = primal.states(7,:);


% alphadot  = primal.states(8,:);

% alphadot2  = primal.controls(1,:)*scale.gammadot; %

% alphadot  = primal.controls(1,:)*scale.gammadot; %
% etadot  = primal.controls(2,:);

alphadot = 0*ones(1,length(alpha));
etadot  = primal.controls(1,:);

time = primal.nodes;
% %=========================================================================================



global q
global M
global Fd
global L
[Vdot,xidot,phidot,gammadot,a,zetadot, q, M, Fd, rho,L] = VehicleModelReturn(time, gamma, V, v, nodes,scattered, Atmosphere,zeta,phi,xi,alpha,eta);

%======================================================

% XDOT = [Vdot;a; gammadot; alphadot; zetadot; phidot; xidot; alphadot2];
% XDOT = [Vdot;a; gammadot; alphadot; zetadot; phidot; xidot; etadot];
XDOT = [Vdot;a; gammadot; zetadot; phidot; xidot; etadot];
end

%======================================================
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

alpha = primal.states(4,:);

zeta = primal.states(5,:);

phi = primal.states(6,:);

xi = primal.states(7,:);

eta = primal.states(8,:);


% alphadot  = primal.states(8,:);

% alphadot2  = primal.controls(1,:)*scale.gammadot; %
alphadot  = primal.controls(1,:)*scale.gammadot; %
etadot  = primal.controls(2,:);
time = primal.nodes;
% %=========================================================================================

global q
global M
global Fd
global L
[Vdot,xidot,phidot,gammadot,a,zetadot, q, M, Fd, rho,L] = VehicleModelReturn(time, gamma, V, v, nodes,scattered, Atmosphere,zeta,phi,xi,alpha,eta);

%======================================================

% XDOT = [Vdot;a; gammadot; alphadot; zetadot; phidot; xidot; alphadot2];
XDOT = [Vdot;a; gammadot; alphadot; zetadot; phidot; xidot; etadot];
end

%======================================================
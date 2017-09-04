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
global interp
global Atmosphere


V = primal.states(1,:)*scale.V ; 

v = primal.states(2,:)*scale.v ; 

vdot = a;
 
gamma = primal.states(3,:)*scale.gamma ; 

alpha = primal.states(4,:);
% alpha = deg2rad(7)*ones(1,length(V));

% zeta = primal.states(5,:);
% 
% phi = primal.states(6,:);
% 
% xi = primal.states(7,:);

% eta = primal.states(8,:);
eta = zeros(1,length(V));


% alphadot  = primal.states(8,:);

% alphadot2  = primal.controls(1,:)*scale.gammadot; %

alphadot  = primal.controls(1,:)*scale.gammadot; %
% etadot  = primal.controls(2,:);

% alphadot = 0*ones(1,length(alpha));
% etadot  = primal.controls(1,:);

time = primal.nodes;
% %=========================================================================================


global q
global M
global Fd
global L
global initial
global phi
[Vdot,xi,phi,gammadot,a,zeta, q, M, Fd, rho,L] = VehicleModelReturn(time, gamma, V, v, nodes,interp, Atmosphere,initial.zeta,initial.phi,initial.xi,alpha,eta);

%======================================================

% XDOT = [Vdot;a; gammadot; alphadot; zetadot; phidot; xidot];

 XDOT = [Vdot;a; gammadot; alphadot];
end

%======================================================
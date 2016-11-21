function XDOT = TwoStage2d(primal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global a
global Fueldt
global const
global scale


V = primal.states(1,:)*scale.V ; 

v = primal.states(2,:)*scale.v ; 

vdot = a;
 
theta = primal.states(3,:)*scale.theta ; 

mfueldot = -Fueldt ; 

thetadot = primal.states(5,:)*scale.thetadot;

% thetadot  = primal.controls(1,:)*scale.theta; %

omegadot  = primal.controls(1,:)*scale.thetadot; %
% %=========================================================================================

Vdot = v.*sin(theta);

%======================================================

XDOT = [Vdot/scale.V;vdot/scale.v; thetadot/scale.thetadot; mfueldot/scale.m; omegadot/scale.thetadot];

end

%======================================================
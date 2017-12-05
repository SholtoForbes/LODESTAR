function XDOT = ThirdStageDynamics(primal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


V = primal.states(1,:); 

v = primal.states(2,:) ; 
 
gamma = primal.states(3,:) ; 

alpha = primal.states(4,:);


alphadot  = primal.controls(1,:); %


time = primal.nodes;
% %=========================================================================================



global rdot
global vdot
global gammadot
global mdot
%======================================================

rdot(end+1) = 0;
vdot(end+1) = 0;
gammadot(end+1) = 0;

XDOT = [rdot;vdot; gammadot; alphadot];

end

%======================================================
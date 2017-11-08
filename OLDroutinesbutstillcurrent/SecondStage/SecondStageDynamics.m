function XDOT = TwoStage2d(primal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global a
global zetadot
global Engine
global Vehicle

v = primal.states(2,:); 
vdot = a;
gamma = primal.states(3,:); 
mfueldot = -Engine.Fueldt; 
gammadot = primal.states(5,:);
omegadot  = primal.controls(1,:); 

%==========================================================================

Vdot = v.*sin(gamma);

%==========================================================================

XDOT = [Vdot;vdot; gammadot; mfueldot; omegadot; zetadot];

end

%======================================================
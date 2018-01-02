function xdot = LanderDynamics(primal)
% Dynamics function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS

h = primal.states(1,:);		
v = primal.states(2,:);		
m = primal.states(3,:);		

Thrust = primal.controls;

% Obviously all this reassignment is not necessary but it's 
% easy to debug the code.  Some penalty in computational 
% performance is incurred

%======================================================================
% Equations of Motion:
%======================================================================
hdot = v;
vdot = - CONSTANTS.g + Thrust./m;
mdot = - Thrust./CONSTANTS.ve;
%=======================================================================
% creating hdot, vdot, mdot is inefficient computing but easier to debug
%=======================================================================
xdot = [hdot; vdot; mdot];

% all done! 
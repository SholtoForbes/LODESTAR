function XDOT = BadBrac1Dynamics(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics for the Bad Brachistochrone Problem
% Template for A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global CONSTANTS UNITS SCALES

xBar = primal.states(1,:);	    
yBar = primal.states(2,:);		
vBar = primal.states(3,:);		

thetaBar  = primal.controls;

%==========================================================
% Equations of Motion:
%==========================================================
xdotBar = SCALES.xdot*vBar.*sin(thetaBar*UNITS.theta);
ydotBar = SCALES.ydot*vBar.*cos(thetaBar*UNITS.theta);							 
vdotBar = SCALES.vdot*CONSTANTS.g*cos(thetaBar*UNITS.theta);
%===========================================================
XDOT = [xdotBar; ydotBar; vdotBar];
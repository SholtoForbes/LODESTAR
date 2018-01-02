function XDOT = Brac1Dynamics(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics for the Brac:1 Formulation of the Brachistochrone Prob
% Template for A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = primal.states(1,:);	    
y = primal.states(2,:);		
v = primal.states(3,:);		


theta  = primal.controls;

%=======================================================
% Equations of Motion:
%=======================================================
xdot = v.*sin(theta);
ydot = v.*cos(theta);							 
vdot = 9.8*cos(theta);  
%======================================================
XDOT = [xdot; ydot; vdot];
%======================================================
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
gamma = primal.states(4,:);	
alpha = primal.states(5,:);	
x = [h; v ;m ;gamma; alpha];

u = primal.controls;
phase = 'postpitch';
global scattered


dz = rocketDynamics(x,u,phase,scattered);

xdot = dz;

% all done! 
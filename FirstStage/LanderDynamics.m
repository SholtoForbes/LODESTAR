function xdot = LanderDynamics(primal)
% Dynamics function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS
global AOAScale
h = primal.states(1,:);		
v = primal.states(2,:);		
m = primal.states(3,:);	
gamma = primal.states(4,:);	
alpha = primal.states(5,:)/AOAScale;	
x = [h; v ;m ;gamma; alpha];

u = primal.controls/AOAScale;
phase = 'postpitch';
global scattered

global q
[dz,q] = rocketDynamics(x,u,phase,scattered);
dz(5) = dz(5)*AOAScale;
xdot = dz;

% all done! 
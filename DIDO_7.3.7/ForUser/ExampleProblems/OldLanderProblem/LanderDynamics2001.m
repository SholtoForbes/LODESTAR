function residuals = LanderDynamics2001(primal)
% Dynamics function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS

h = primal.states(1,:);		hdot = primal.statedots(1,:);
v = primal.states(2,:);		vdot = primal.statedots(2,:);
m = primal.states(3,:);		mdot = primal.statedots(3,:);

Thrust = primal.controls;

% Obviously all this reassignment is not necessary but it's 
% easy to debug the code.  Some penalty in computational 
% performance is incurred

% preallocate residuals for good MATLAB computing.  This can 
% be done by writing:
%
% [Nx, Nn]  = size(primal.states);
% residuals = zeros(Nx, Nn);
%
% A better idea due to Prof. C. D. Hall is to obviate the need to 
% to repeatedly compute "size" and "zeros" by simply allocating:

residuals = primal.states;


%=======================================================
% Equations of Motion:
residuals(1,:) = hdot - v;
residuals(2,:) = vdot + CONSTANTS.g - Thrust./m;
residuals(3,:) = mdot + Thrust./CONSTANTS.ve;
%=======================================================

% Note: A neat trick for preallocation is to write the last equation first.
% MATLAB will automatically preallocate the remainder of the rows!

% all done! 
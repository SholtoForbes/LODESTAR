function residuals = stageDynamics(primal)
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global CONSTANTS

r = primal.states(1,:);   rdot = primal.statedots(1,:);
v = primal.states(2,:);   vdot = primal.statedots(2,:);
m = primal.states(3,:);   mdot = primal.statedots(3,:);

Thrust = primal.controls;

[Nx, Nn] = size(primal.states);   residuals = zeros(Nx, Nn);


%  Equations of motion:
%============================================================
%
residuals(1,:) = rdot - v;
residuals(2,:) = vdot - Thrust./m  + 1./r.^2;
residuals(3,:) = mdot + Thrust./CONSTANTS.Ve;
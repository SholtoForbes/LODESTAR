function xdot = stageDynamics(primal)
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Example 2, pp.403-404 in I.M.Ross and F. Fahroo,
% "Pseudospectral Knotting Methods for Solving Optimal Control Problems,"
% Journal of Guidance, Control and Dynamics, Vol.27, No.3, May-June 2004.
%-----------------------------------------------%
global CONSTANTS

r = primal.states(1,:);   
v = primal.states(2,:);   
m = primal.states(3,:); 

Thrust = primal.controls;

%============================================================
%  Equations of motion:
%============================================================
%
rdot =  v;
vdot = Thrust./m  - 1./r.^2;
mdot = -Thrust./CONSTANTS.Ve;

xdot = [rdot; vdot; mdot];
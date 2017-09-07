%------------------------- Bryson-Denham Problem -------------------------%
% This problem is taken from the following reference:                     %
% Bryson, A. E. and Ho, Y-C, "Applied Optimal Control:  Optimization,     %
% Estimation, and Control," Hemisphere Publishing, 1975.                  %
%-------------------------------------------------------------------------%
clear all; clc

%--------------------------------------------------------------------------%
%---------------- Set Up Auxiliary Data for Problem -----------------------%
%--------------------------------------------------------------------------%
auxdata.T  = 0.1405;
auxdata.m0 = 1;
auxdata.dm = 0.0749;
auxdata.mu = 1;

%--------------------------------------------------------------------------%
%--------------- Set Up Bounds on State, Control, and Time ----------------%
%--------------------------------------------------------------------------%
umax        = 0.002;
t0          = 0;
tfMin       = 0;
tfMax       = 500;
r0          = 1; 
theta0      = 0; 
vr0         = 0;
vtheta0     = 1; 
rf          = 4;
vrf         = 0;
vthetaf     = 0.5;
rmin        = 1; 
rmax        = 10;
thetamin    = 0; 
thetamax    = 200*pi;
vrmin       = -100;
vrmax       = +100;
vthetamin   = -100; 
vthetamax   = +100;
u1min       = -umax;
u1max       = +umax;
u2min       = -umax; 
u2max       = +umax;
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [r0, theta0, vr0, vtheta0];
bounds.phase.initialstate.upper = [r0, theta0, vr0, vtheta0];
bounds.phase.state.lower = [rmin, thetamin, vrmin, vthetamin];
bounds.phase.state.upper = [rmax, thetamax, vrmax, vthetamax];
bounds.phase.finalstate.lower = [rf, thetamin, vrf, vthetaf];
bounds.phase.finalstate.upper = [rf, thetamax, vrf, vthetaf];
bounds.phase.control.lower = [u1min, u2min];
bounds.phase.control.upper = [u1max, u2max];
% bounds.phase.path.lower  = 1;
% bounds.phase.path.upper  = 1;
% bounds.eventgroup.lower = [0];
% bounds.eventgroup.upper = [0];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Guess ---------------------------%
%--------------------------------------------------------------------------%
tGuess      = [t0; tfMax];
rGuess      = [r0; rf];
thetaGuess  = [theta0; 20*pi];
vrGuess     = [vr0; vrf];
vthetaGuess = [vtheta0; vthetaf];
u1Guess   = [0;0];
u2Guess   = [0;0];
guess.phase.time    = [tGuess];
guess.phase.state   = [rGuess, thetaGuess, vrGuess, vthetaGuess];
guess.phase.control = [u1Guess, u2Guess];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Mesh ----------------------------%
%--------------------------------------------------------------------------%
N = 10;
meshphase.colpoints = 4*ones(1,N);
meshphase.fraction   = ones(1,N)/N;

%--------------------------------------------------------------------------%
%-------------------------- Set Up for Solver -----------------------------%
%--------------------------------------------------------------------------%
setup.name = 'Ross-Low-Thrust-Problem';
setup.functions.continuous = @lowThrustRossContinuous;
setup.functions.endpoint = @lowThrustRossEndpoint;
setup.displaylevel = 2;
setup.nlp.solver = 'ipopt';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.mesh.method = 'hp-LiuRao-Legendre';
setup.mesh.tolerance = 1e-6;
setup.mesh.phase = meshphase;
setup.scales.method = 'none';

%--------------------------------------------------------------------------%
%-------------------- Solve Problem and Extract Solution ------------------%
%--------------------------------------------------------------------------%
output = gpops2(setup);

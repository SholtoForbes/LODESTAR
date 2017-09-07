%---------------------- Double-Integrator Problem ------------------------%
% This example can be found in any good book on optimal control theory    %
%-------------------------------------------------------------------------%
clear all; close all; clc

auxdata.rho = 2.5;
auxdata.k   = 1.5;
%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0Min = 0;  t0Max = 0;
tfMin = 2;  tfMax = 2;
x10   = +1; x1Min = -10; x1Max = -x1Min;
x20   = +0.01; x2Min = -10; x2Max = -x2Min;
uMin = +0.1;   uMax = 0.5; % (0.4, 0.3)

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0Min;
bounds.phase.initialtime.upper = t0Max;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [x10, x20];
bounds.phase.initialstate.upper = [x10, x20];
bounds.phase.state.lower = [x1Min, x2Min];
bounds.phase.state.upper = [x1Max, x2Max];
bounds.phase.finalstate.lower = [x1Min, x2Min];
bounds.phase.finalstate.upper = [x1Max, x2Max];
bounds.phase.control.lower = [uMin];
bounds.phase.control.upper = [uMax];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [t0Min; tfMax];
x1Guess              = [x10; x10];
x2Guess              = [x20; x20];
uGuess               = [uMax; uMax];
guess.phase.state    = [x1Guess, x2Guess];
guess.phase.control  = [uGuess];
guess.phase.time     = [tGuess];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao';
mesh.tolerance       = 1e-6;
mesh.maxiterations    = 20;
mesh.colpointsmin    = 4;
mesh.colpointsmax    = 10;
mesh.phase.colpoints = 4;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.mesh                        = mesh;
setup.name                        = 'Chemical-Reaction';
setup.auxdata                     = auxdata;
setup.functions.endpoint          = @chemicalReactionEndpoint;
setup.functions.continuous        = @chemicalReactionContinuous;
setup.displaylevel                = 2;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%----------------------- Solve Problem Using GPOPS2 ----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);

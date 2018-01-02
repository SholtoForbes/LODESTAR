%  StageMain:  Main File (Script) for the Staging Problem
%======================================================================
% Example file for DIDO 7.3
% TBD DIDO User's Manual
% Jim B. Ross & I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;						% very good idea to begin with this
%==============================================================
global CONSTANTS

CONSTANTS.Tmax     = 3.5;      % Normalized Ultra Max thrust
CONSTANTS.Ve       = 0.5;      % Normalized Isp*g0

%========================================================================
%  Define the problem:
%========================================================================

stageProb.cost      = 'stageCost';
stageProb.dynamics  = 'stageDynamics2001';
stageProb.events    = 'stageEvents';
%----- comment and uncomment the following -----
stageProb.path      = 'stagePath';
%----------------------------------------------
%========================================================================
%  Setup the knots and some more constants
%========================================================================

t0 = 0;   t1Guess = .04;   tfMax = 100;   tfGuess = .2;

t2Guess =0.155

knots.locations    = [t0  t1Guess t2Guess tfGuess];
knots.definitions  = {'hard','hard','soft', 'hard'};
knots.bounds.lower = [0 eps 0.15 .06];           
knots.bounds.upper = [0 tfMax/2 0.25 tfMax];
knots.numNodes     = [10 24 6];

Ntotal = sum(knots.numNodes);

%========================================================================
%  Set up the bounds and boundary conditions
%========================================================================

r0 = 1.0;    v0 = 0;     m0 = 1.0;     % initial conditions
dr1 = 0;    dv1 = 0;    stg1Fuel = 0.2;  stg2m0 = .7;      % event conditions  
mf = 0.5;                              % boundary condition


bounds.lower.states = [1; -20; mf];
bounds.upper.states = [20; 20; m0];

%----------Segmented Controls Example -----------------------
bounds.lower.controls = [0, 0];;
bounds.upper.controls = [CONSTANTS.Tmax/2 CONSTANTS.Tmax/3];
%--------- values arbitrarily chosen for illustration --------

bounds.lower.events = [r0; v0; m0; dr1; dv1; 0; stg2m0; mf];
bounds.upper.events = [r0; v0; m0; dr1; dv1; stg1Fuel; stg2m0; mf];

%-------Segmented Bounds Example ----------------------------
bounds.lower.path = [0 0];
bounds.upper.path = [CONSTANTS.Tmax/3 CONSTANTS.Tmax/4]; % arbitrarily chosen
%-------Comment and uncomment this with "path function" -----


%========================================================================
%  Provide a guess
%========================================================================

rfGuess = 1.3;

% Initial guess values
%========================================================================
guess.states(1,:) = [r0, 1.15, rfGuess];
guess.states(2,:) = [v0, .1,  0];
guess.states(3,:) = [m0, .6,  mf];
guess.controls    = [3.5, 3.5, 0];
guess.time        = [t0, .125, tfGuess];

%% call dido

%***************************************
% NEW COMPATIBILITY REQUIREMENT
%***************************************
stageProb.format       = '2001';
%***************************************

[cost, primal] = dido(stageProb, knots, bounds, guess);

save stageOutSoft

%============================================================================
figure;
plot(primal.nodes, primal.states, '*', primal.nodes, primal.controls, '+');
legend('altitude', 'speed', 'mass', 'thrust');
xlabel('normalized time units');
ylabel('normalized units');
hold on;
plot(primal.nodes, primal.states, primal.nodes, primal.controls);
hold off;
%=============================================================================
%  StageProblem:  Problem or Main File (Script) for the Staging Problem
%======================================================================
% Example file for DIDO
% TBD DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Example 2, pp.403-404 in I.M.Ross and F. Fahroo,
% "Pseudospectral Knotting Methods for Solving Optimal Control Problems,"
% Journal of Guidance, Control and Dynamics, Vol.27, No.3, May-June 2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------
clear all;						% very good idea to begin with this
%==================================================================

global CONSTANTS

CONSTANTS.Tmax     = 3.5;      % Normalized Ultra Max thrust
CONSTANTS.Ve       = 0.5;      % Normalized Isp*g0

%===================
% Problem variables:
%-------------------
% states = (r, v, m)
% controls = Thrust
%===================

%---------------------------------------
% bounds the state variables
%---------------------------------------

r0 = 1.0;    v0 = 0;     m0 = 1.0;     % initial conditions
mf = 0.5;                              % boundary condition


bounds.lower.states = [1; -20; mf];
bounds.upper.states = [20; 20; m0];

%-------------------------------------------------------------
% bounds the control variables: illustrating segmented bounding
%-------------------------------------------------------------

bounds.lower.controls = [0, 0];;
bounds.upper.controls = [CONSTANTS.Tmax/2 CONSTANTS.Tmax/3];

% The above bounds are actually non-binding constraints; they may be set to
% binding constraints but written here as non-binding to illustrate how to
% do segmented bounding. The binding constraints are in the path funtion for
% this problem; again, for the purposes of illustration.


%-------------------------
% bound the time intervals
%-------------------------

t0 = 0;   tfMax = 100;  

% bound the interior event time between eps and half of maxtf; the 0.06 is
% a swag.

bounds.lower.time = [0 eps .06];           
bounds.upper.time = [0 tfMax/2 tfMax];  


%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------
% See events file for definition of events function and Eqs 63-67 in the
% Ross-Fahroo paper

dr1 = 0;    dv1 = 0;    stg1Fuel = 0.2;  stg2m0 = .7;      % event conditions data 

bounds.lower.events = [r0; v0; m0; dr1; dv1; 0; stg2m0; mf];
bounds.upper.events = [r0; v0; m0; dr1; dv1; stg1Fuel; stg2m0; mf];


%-------Illustrating Segmented Bounds  ----------------------------
% bounds.lower.path = [0 0];
% bounds.upper.path = [CONSTANTS.Tmax/3 CONSTANTS.Tmax/4]; 

bounds.lower.path = [0 0];
bounds.upper.path = [0.07 0.07]; 


% See data in the Ross-Fahroo paper
%-------This is binding.  If these values are put in the box constraints
% forf the controls, then a path function is not necessary; this is done 
% merely to illustrate how to do segmented bounding -----


%========================================================================
%  Define the problem using DIDO Expressions:
%========================================================================

stageProb.cost      = 'stageCost';
stageProb.dynamics  = 'stageDynamics';
stageProb.events    = 'stageEvents';
%----- comment and uncomment the following if box constraints on the thrust is binding -----
stageProb.path      = 'stagePath';
%----------------------------------------------

stageProb.bounds = bounds;

%========================================================================
%  Tell DIDO that all your events are 'hard'; this redundant statement will
%  be removed in future versions of DIDO
%========================================================================
algorithm.knots.definitions  = {'hard','hard','hard'};
%========================================================================
% For a 3 stage problem, you need to add another 'hard' etc. Yes, it's
% redundant.

%=======================================================================
% Tell DIDO how many nodes you want to use for the first and second stages
%=======================================================================
algorithm.nodes     = [10 30];

%======================================================================
% Providea Guess.  The Guess-Free Option is not available for "Knots"
%====================================================================

t1Guess = .04;  % Guess for the interior event time
tfGuess = .2;   % Guess for the final time

% Tell DIDO these guesses using the "knots" field in algorithm:

algorithm.knots.locations    = [t0  t1Guess tfGuess];

% State and control guesses need to include a guess for the interior
% variables

rfGuess = 1.3;

%========================================================================
guess.states(1,:) = [r0, 1.15, rfGuess];
guess.states(2,:) = [v0, .1,  0];
guess.states(3,:) = [m0, .6,  mf];
guess.controls    = [3.5, 3.5, 0];
guess.time        = [t0, .125, tfGuess];
%=======================================================================

% Tell DIDO the guess.  Note: The guess-free option is not available when
% using "knots"
%========================
algorithm.guess = guess;
%========================

%% call dido

[cost, primal] = dido(stageProb, algorithm);

save stageOut

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


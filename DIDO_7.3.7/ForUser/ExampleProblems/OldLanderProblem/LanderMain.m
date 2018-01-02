% LanderMain: Main File (Script) for the Moon-Landing Problem
%==============================================================
% Straightforward Setup
%--------------------------------------------------------------
% Example file for DIDO in 2001 format
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;					% very good idea to begin with this
%==============================================================
global CONSTANTS

CONSTANTS.g     = 1.000;     % g of moon
CONSTANTS.ve    = 2.349;     % exhaust velocity (very poor man's rocket!)
CONSTANTS.Tmax  = 1.227;     % all are in normalized units

%-------------------------
% Define the problem:
%-------------------------
moonLanderStraight.cost 		= 'LanderCost';
moonLanderStraight.dynamics	    = 'LanderDynamics2001';
moonLanderStraight.events		= 'LanderEvents';		
% no path constraints

%***************************************
% NEW COMPATIBILITY REQUIREMENT
%***************************************
moonLanderStraight.format       = '2001';
%***************************************

%----------------------------------------
% Setup the knots and some more constants
%----------------------------------------

t0			= 0;
tfMax 	    = 20;     % large upper bound; better than choosing Inf
tfGuess 	= 2;			
% arbitrary guess of final time (for the scaled problem!)

knots.locations		= [t0 tfGuess];		% time free problem
knots.definitions	= {'hard', 'hard'};
knots.bounds.lower 	= [0 0];				
knots.bounds.upper	= [0 tfMax];			 
knots.numNodes		= [30];					  

% Note the fixed hard knot at t0 and a free hard knot at tf

Nn	= sum(knots.numNodes);

% total number of nodes; represents accuracy over the time span

%--------------------------------------------------------------
% Set up the bounds and boundary conditions in the SCALED units
%--------------------------------------------------------------


h0 = 1.000;     v0 = -0.783;    m0 = 1.000;     % initial conditions
hf = 0.0;       vf =  0.0;      mfmin = 0.01;   % boundary conditions

% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.

bounds.lower.states = [-20; -20; mfmin];
bounds.upper.states = [ 20;  20; m0];

bounds.lower.controls = 0;
bounds.upper.controls = CONSTANTS.Tmax;

bounds.lower.events = [h0; v0; m0; hf; vf];	
bounds.upper.events = bounds.lower.events;

%-----------------------------
% Provide a guess
%-----------------------------
mfguess = 0.1; 

guess.states(1,:)	= [h0, hf];
guess.states(2,:)	= [v0, vf];
guess.states(3,:)	= [m0, mfguess];
guess.controls		= [0.01, 0.01];
guess.time			= [t0, tfGuess];

% DIDO will interpret this to be a straight line passing through the 
% specified points. Normally you should specify a "better" guess using 
% as many points as you can.  See the DIDO User's Manual on generating good guesses.

startTime = cputime;
[cost, primal, dual] = dido(moonLanderStraight, knots, bounds, guess);
runTime = cputime-startTime

save landerOut


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
% LanderMainSoft: Main File for the Moon-Landing Problem
%============================================================
% With Soft Knots
%------------------------------------------------------------
% Example file for DIDO 
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;					% very good idea to begin with this
%==============================================================
global CONSTANTS

CONSTANTS.g     = 1.000;     % g of moon
CONSTANTS.ve    = 2.349;     % exhaust velocity (very poor man's rocket!)
CONSTANTS.Tmax  = 1.227;     % all are in normalized units

%--------------------
% Define the problem
%--------------------

moonLanderSoft.cost 	= 'LanderCost';
moonLanderSoft.dynamics	= 'LanderDynamics2001';

%***************************************
% NEW COMPATIBILITY REQUIREMENT
%***************************************
moonLanderSoft.format       = '2001';
%***************************************
%------
% Knots
%------

t0			= 0;
tfMax 	    = 1000;
tfGuess 	= 1.5;			
% arbitrary guess of final time (normalized units)

knots.locations		= [t0 0.3 tfGuess];
knots.definitions	= {'hard', 'soft' 'hard'};
knots.bounds.lower 	= [0 0.1 0];
knots.bounds.upper	= [0 0.6 tfMax];
knots.numNodes		= [4 12];

% Look: same number of total points spread unevenly.  
% Try ruNning this with an even spread of points.  Explain what you see.

Nn	= sum(knots.numNodes);								
% total number of nodes; represents accuracy over the time span

%---------------------------------------
% BOUNDS: including "Events" and "Path"
%---------------------------------------

h0 = 1.000;     v0 = -0.783;    m0 = 1.000;     % initial conditions
hf = 0.0;       vf =  0.0;      mfmin = 0.01;   % boundary conditions


bounds.lower.states(1,:) = [h0, -20,  hf];	
bounds.upper.states(1,:) = [h0,  20,  hf];	 
bounds.lower.states(2,:) = [v0, -20,  vf];	   
bounds.upper.states(2,:) = [v0,  20,  vf];	
bounds.lower.states(3,:) = [m0, mfmin, mfmin];	
bounds.upper.states(3,:) = [m0, m0, m0];

% These are *NOT* the minimum requirements for specifying bounds.
% The "middle" point (e.g. -20) may be omitted, in this case,
% DIDO will assume unbounded interior points.

bounds.lower.controls = 0;
bounds.upper.controls = CONSTANTS.Tmax;


%------------------------------------
% load a "good guess" for soft knots 
%------------------------------------

load landerOut primal
guess.states 	= primal.states;
guess.controls	= primal.controls;
guess.time		= primal.nodes;
knots.locations(1) = primal.nodes(1); knots.locations(end) = primal.nodes(end);

%------
% GUESS
%------
%guess.states(1,:)	= [h0, hf];				 
%guess.states(2,:)	= [v0, vf];				 
%guess.states(3,:)	= [m0, mfguess];		  
%guess.controls		= [0.01, 0.01];		
%guess.time			= [t0, tfGuess];		



[cost, primal, dual] = dido(moonLanderSoft, knots, bounds, guess);

save LanderOutSoft

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

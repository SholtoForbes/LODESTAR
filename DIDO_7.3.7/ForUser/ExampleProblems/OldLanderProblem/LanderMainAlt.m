% LanderMainAlt  Main File (Script) for the 1-DOF Moon-Landing Problem 
% Alternative Setup: No Events or Path file used
% ==================================================================
% Example file for DIDO  
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;				% we recommend you always do this

global CONSTANTS

CONSTANTS.g 	= 1.0;		% in normalized units
CONSTANTS.Isp	= 1.0;		% ditto
CONSTANTS.Tmax	= 0.5;		% ditto
CONSTANTS.ve	= 1.0;		% ditto

%--------------------
% Define the problem
%--------------------
%***************************************
% NEW COMPATIBILITY REQUIREMENT
%***************************************
moonLander.format       = '2001';
%***************************************

moonLander.cost 	= 'LanderCostMayer';
moonLander.dynamics	= 'LanderDynamics2001';


% Events and path specified "indirectly" in bounds

%--------
%% KNOTS:
%--------

t0			= 0;
tfMax 	= Inf;
tfGuess 	= 1.5;			
% arbitrary guess of final time (normalized units)

knots.locations		= [t0 tfGuess];
knots.definitions	= {'hard', 'hard'};
knots.bounds.lower 	= [0 0];
knots.bounds.upper	= [0 tfMax];
knots.numNodes		= [32];

Ntotal	= sum(knots.numNodes);								
% total number of nodes; represents accuracy over the time span

%---------------------------------------
% BOUNDS: including "Events" and "Path"
%---------------------------------------

h0=1.0; v0=-0.05; m0=1.0;  		% initial conditions
hf=0.0; vf= 0.0;   					% final conditions
mfguess = 0.1; mfmin = 0.01;		% guess of missing boundary condition

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

%------
% GUESS
%------
guess.states(1,:)	= [h0, hf];				 
guess.states(2,:)	= [v0, vf];				 
guess.states(3,:)	= [m0, mfguess];		  
guess.controls		= [0.01, 0.01];		
guess.time			= [t0, tfGuess];		


[cost, primal, dual] = dido(moonLander, knots, bounds, guess);

lam_v = dual.dynamics(2,:);
lam_m = dual.dynamics(3,:);
m = primal.states(3,:);

save LanderOutAlt

%============================================================================
figure;
plot(primal.nodes, primal.states, '*', primal.nodes, primal.controls, '+');
legend('altitude', 'speed', 'mass', 'thrust');
xlabel('normalized time units');
ylabel('normalized units');
hold on;
plot(primal.nodes, primal.states, primal.nodes, primal.controls);
hold off;

% plot the switching function
figure;
S_M = lam_v./m - lam_m/CONSTANTS.ve;
S_L = 1 + S_M; % lagrange cost switch function
plot(primal.nodes, S_M, primal.nodes, primal.controls,'+');
%=============================================================================
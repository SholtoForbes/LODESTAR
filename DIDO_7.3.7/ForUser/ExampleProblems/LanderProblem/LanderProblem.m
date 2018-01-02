% LanderProblem: Problem File (Script) for the Moon-Landing Problem
%==============================================================
% Minimalist's Setup
%--------------------------------------------------------------
% Example file for DIDO 
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;					% very good idea to begin with this
%==============================================================
global CONSTANTS

CONSTANTS.g     = 1.000;     % g of moon
CONSTANTS.ve    = 2.349;     % exhaust velocity (very poor man's rocket!)
CONSTANTS.Tmax  = 10;     % all are in normalized units

%-----------------------------------
% Define the problem function files:
%-----------------------------------
MoonLander.cost 		= 'LanderCost';
MoonLander.dynamics	    = 'LanderDynamics';
MoonLander.events		= 'LanderEvents';	
% MoonLander.Path		= 'LanderPath';
% no path constraints

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 1.8;     % large upper bound; do not choose Inf
% arbitrary guess of final time (for the scaled problem!)

bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 

h0 = 0.0;     v0 = 0;    m0 = 1.000;     % initial conditions
hf = 1.0;       vf =  1.0;      mfmin = 0.01;   % boundary conditions

% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.

bounds.lower.states = [-20; -20; mfmin];
bounds.upper.states = [ 20;  20; m0];

bounds.lower.controls = 0;
bounds.upper.controls = CONSTANTS.Tmax;

bounds.lower.events = [h0; v0; m0; hf];	
bounds.upper.events = bounds.lower.events;

bounds.lower.Path = -0.8;
bounds.upper.Path = 1.0;


%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

MoonLander.bounds = bounds;

%------------------------------------------------------
% Select the number of nodes for the spectral algorithm
%------------------------------------------------------

algorithm.nodes = [30];  % somewhat arbitrary number; theoretically, the 
                         % larger the number of nodes, the more accurate 
                         % the solution (but, practically, this is not
                         % always true!)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At this point you can now run dido several different ways %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====================================
% Recommended approachs for beginners
%====================================
%------------------------------------------------------------------
% Exercise # 1: Compare the run times for different calling options
%------------------------------------------------------------------

approach = 1

switch approach
    case 1

%---------------------------
% call DIDO: Simplest option
%---------------------------

startTimeNoGuess = cputime;
[cost, primal, dual] = dido(MoonLander, algorithm);
noGuessRunTime = cputime - startTimeNoGuess

%--- save data ---
save NoGuessLanderOutput


%--- plot data ---

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

%============================================================================
figure;
plot(primal.nodes, dual.dynamics, '-*');
legend('\lambda_h', '\lambda_v', '\lambda_m');
xlabel('normalized time units');
ylabel('normalized units');

%=============================================================================
return

    case 2

%---------------------------
% call DIDO: Guess option
%---------------------------      

% Use your knowledge of the problem to give DIDO some starting point

mfguess     = 0.1; 
t0			= 0;
tfGuess 	= 1.5;			
% slightly educated guess of final time (for the scaled problem!)

guess.states(1,:)	= [h0, hf];
guess.states(2,:)	= [v0, vf];
guess.states(3,:)	= [m0, mfguess];
guess.controls		= [0.01, 0.01];
guess.time			= [t0, tfGuess];

% This is a "two point" guess. Normally you should specify a "better" guess using 
% as many points as you can.  See the DIDO User's Manual for generating good guesses.

%----------------------------
% Tell DIDO you have a guess:

algorithm.guess = guess;
%----------------------------

startTimeUserguess = cputime;
[cost, primal, dual] = dido(MoonLander, algorithm);
runTimeWguess = cputime-startTimeUserguess


save UserGuessLanderOut


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



end
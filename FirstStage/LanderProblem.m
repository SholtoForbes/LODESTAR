


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
clc
%==============================================================


%-----------------------------------
% Define the problem function files:
%-----------------------------------
MoonLander.cost 		= 'LanderCost';
MoonLander.dynamics	    = 'LanderDynamics';
MoonLander.events		= 'LanderEvents';	

global scattered

addpath TrajOpt-master

Aero = dlmread('FirstStageAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

global SPARTANscale
% SPARTANscale = 0.75
SPARTANscale = 1

% mRocket = 27000; %(kg)  %Total lift-off mass
mRocket = 27100; %(kg)  %Total lift-off mass
mFuel = 0.89*mRocket;  %(kg)  %mass of the fuel

mSpartan = 8755.1*SPARTANscale;

mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = 0;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan = [0 15];
y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0];
% [t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,Tmax,phase), tspan, y0);
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,phase,scattered), tspan, y0);  

% FOR TESTINGm, see where it gets 
phase = 'postpitch';
Tratio = 1;
tspan = [0 mFuel/156]; % flies for way too long
postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) deg2rad(0)];
[t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,0,phase,scattered), tspan, postpitch0);

y
postpitch
postpitch(end,4)
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
global AOAScale
AOAScale = 1;
h0 = y(end,1);  
v0 = y(end,2);  %
m0 = y(end,3);  
gamma0 = deg2rad(89.9);    % pitchover 
% gamma0 = deg2rad(89.999);    % pitchover 

vF = 1850;  
mF = mEmpty+mSpartan;  %Assume that we use all of the fuel
gammaF = deg2rad(1);
% hF = 26550;
% hF = 45000;
hF = 26000;

hLow = 0;   %Cannot go through the earth
hUpp = 70000;  

vLow = 0; 
vUpp = 3000;  

mLow = mEmpty;
mUpp = mTotal;

gammaLow = deg2rad(-1);
gammaUpp = deg2rad(89.9);



uLow = [-.005]*AOAScale; % Can do either AoA or thrust
uUpp = [.005]*AOAScale; 

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 300;     % large upper bound; do not choose Inf


bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 


% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.


bounds.lower.states = [hLow; vLow; mF-1;gammaLow;-deg2rad(3)*AOAScale];
bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(3)*AOAScale];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;

% bounds.lower.events = [h0; v0; m0; gamma0; hF; mF; gammaF];	
bounds.lower.events = [h0; v0; m0; gamma0; hF; mF];	
% bounds.lower.events = [h0; v0; m0; gamma0; mF; gammaF];	
bounds.upper.events = bounds.lower.events;



%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

MoonLander.bounds = bounds;

%------------------------------------------------------
% Select the number of nodes for the spectral algorithm
%------------------------------------------------------

algorithm.nodes = [50];  % somewhat arbitrary number; theoretically, the 
                         % larger the number of nodes, the more accurate 
                         % the solution (but, practically, this is not
                         % always true!)



    



t0			= 0;
tfGuess 	= 120;			
% slightly educated guess of final time (for the scaled problem!)

guess.states(1,:)	= [h0, hF];
guess.states(2,:)	= [v0, vF];
guess.states(3,:)	= [m0, mF];
guess.states(4,:)	= [gamma0,gammaF];
guess.states(5,:)	= [-deg2rad(0.03)*AOAScale, -deg2rad(0.03)*AOAScale];
guess.controls		= [0.0, 0.0];
guess.time			= [t0, tfGuess];



%----------------------------
% Tell DIDO you have a guess:

algorithm.guess = guess;
%----------------------------

startTimeUserguess = cputime;
[cost, primal, dual] = dido(MoonLander, algorithm);
runTimeWguess = cputime-startTimeUserguess


t = primal.nodes;


%============================================================================
V = primal.states(1,:);
v = primal.states(2,:);
m = primal.states(3,:);
gamma = primal.states(4,:);
alpha = primal.states(5,:);
figure;
hold on
plot(primal.nodes, gamma);
plot(primal.nodes, v/1000);
plot(primal.nodes, V/10000);
plot(primal.nodes, alpha*10);

mu_1 = dual.states(1,:);
mu_2 = dual.states(2,:);
mu_3 = dual.states(3,:);
mu_4 = dual.states(4,:);
mu_5 = dual.states(5,:);

mu_u = dual.controls; % NOTE: This deviates from 0, as the controls are set as a buffer. Do not set a parameter directly tied to the vehicle model as the control.

%GRADIENT NORMALITY CONDITION

% Lagrangian of the Hamiltonian 
dLHdu = dual.dynamics(3,:) + mu_u; % 

figure(5)

plot(t,dLHdu,t,mu_1,t,mu_2,t,mu_3,t,mu_4,t,mu_5,t,mu_u);
legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
title('Validation')
%=============================================================================


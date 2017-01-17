


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
clc
%==============================================================
% global zeta
global phi
global q


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

% mRocket = 27100; %(kg)  %Total lift-off mass
% mRocket = 21000; %(kg)  %Total lift-off mass (this is almost exactly half the mass of a falon 1 first stage, and would give a length of 8.1m scaled to exactly half size (9m with margin of error))
 mRocket = 20000; %(kg)  %Total lift-off mass (this is almost exactly half the mass of a falon 1 first stage, and would give a length of 8.1m scaled to exactly half size (9m with margin of error))
% mFuel = 0.89*mRocket;  %(kg)  %mass of the fuel
% mFuel = 0.939*mRocket;  %(kg)  %mass of the fuel ( from falcon 1 users manual)

mEngine = 470; % Mass of Merlin 1C
% mFuel = 0.939*mRocket; 
mFuel = 0.946*(mRocket-mEngine);  %smf without engine
% mFuel = 0.939*mRocket -mEngine; 
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
tspan = [0 30];
y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0];
% [t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,Tmax,phase), tspan, y0);
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,scattered), tspan, y0);  

% FOR TESTINGm, see where it gets 
phase = 'postpitch';
Tratio = 1;
tspan = [0 mFuel/156]; % flies for way too long
postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) deg2rad(0) 0];
[t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,0,0,phase,scattered), tspan, postpitch0);

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
hF = 26300;

hLow = 0;   %Cannot go through the earth
hUpp = 70000;  

vLow = 0; 
vUpp = 3000;  

mLow = mEmpty;
mUpp = mTotal;

gammaLow = deg2rad(1);
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


bounds.lower.states = [hLow; vLow; mF-1;gammaLow;-deg2rad(3)*AOAScale;0];
bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(3)*AOAScale;2*pi];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;

zetaF = deg2rad(97);

% bounds.lower.events = [h0; v0; m0; gamma0; hF; mF; gammaF];	
% bounds.lower.events = [h0; v0; m0; gamma0; hF; mF; zetaF];	
% bounds.lower.events = [h0; v0; m0; gamma0; mF; gammaF];	

bounds.lower.events = [h0; v0; m0; gamma0; mF; zetaF];	
bounds.upper.events = bounds.lower.events;



%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

MoonLander.bounds = bounds;

%------------------------------------------------------
% Select the number of nodes for the spectral algorithm
%------------------------------------------------------

algorithm.nodes = [70];  % somewhat arbitrary number; theoretically, the 
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
guess.states(6,:)	= [1.35, zetaF];
% guess.states(5,:)	= [0, 0];
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
zeta = primal.states(6,:);
figure;
hold on
plot(primal.nodes, gamma,'color','k','linestyle','-');
plot(primal.nodes, v/1000,'color','k','linestyle','--');
plot(primal.nodes, V/10000,'color','k','linestyle',':');
plot(primal.nodes, rad2deg(alpha)/10,'color','k','linestyle','-.');
legend('Trajectory Angle (degrees/10)','velocity (m/s / 1000)','Altitude (km /10)','angle of attack (rad)')
xlabel('Time (s)')
xlim([0,primal.nodes(end)]);

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

% FOR TESTINGm, see where it gets 
% dalphadt = [diff(alpha)./diff(primal.nodes) 0];
% 
% phase = 'postpitch';
% Tratio = 1;
% tspan = primal.nodes; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) alpha(1)];
% [t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamics(postpitch_f,ControlFunction(t,primal.nodes,dalphadt),phase,scattered), tspan, postpitch0_f);


 phase = 'postpitch';
tspan = primal.nodes; 
postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)];
[t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t,primal.nodes,zeta),ControlFunction(t,primal.nodes,alpha),phase,scattered), tspan, postpitch0_f);


% y
% postpitch_f
% postpitch_f(end,4)

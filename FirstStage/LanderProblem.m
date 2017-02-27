


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
clc
%==============================================================
% global zeta
global phi
global q

global iterative_V
iterative_V = [];
global iterative_t
iterative_t = [];
global iteration
iteration = 1;
global iterative_V_f
iterative_V_f = [];

%-----------------------------------
% Define the problem function files:
%-----------------------------------
MoonLander.cost 		= 'LanderCost';
MoonLander.dynamics	    = 'LanderDynamics';
MoonLander.events		= 'LanderEvents';	
% MoonLander.path		= 'LanderPath';
global scattered

addpath TrajOpt-master

Aero = dlmread('FirstStageAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = scattered.Lift(grid.M,grid.AoA);
grid.Drag = scattered.Drag(grid.M,grid.AoA);
scattered.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
scattered.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');

global SPARTANscale
% SPARTANscale = 0.75
SPARTANscale = 1


mRocket = 18500; % sets the total wet mass of the rocket (first stage only)
 % mFuel = 0.89*mRocket;  %(kg)  %mass of the fuel
% mFuel = 0.939*mRocket;  %(kg)  %mass of the fuel ( from falcon 1 users manual)

mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*mRocket; % Michael said to just use this for simplicity
% mFuel = 0.946*(mRocket-mEngine);  %smf without engine
% mFuel = 0.939*mRocket -mEngine; 
% mSpartan = 8755.1-302+139.4;
mSpartan = 9.9725e+03;

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
tspan = [0 25];
% tspan = [0 40]; 
y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0];
% [t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,Tmax,phase), tspan, y0);
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,scattered), tspan, y0);  

% FOR TESTINGm, see where it gets 
phase = 'postpitch';
Tratio = 1;
tspan = [0 mFuel/156]; % flies for way too long
postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) deg2rad(-1.4) 1.63];
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
% gamma0 = deg2rad(89);    % pitchover 

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



uLow = [-.01]*AOAScale; % Can do either AoA or thrust
uUpp = [.01]*AOAScale; 

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 300;     % large upper bound; do not choose Inf


bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 


% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.


bounds.lower.states = [hLow; vLow; mF-1;-0.1;-deg2rad(4)*AOAScale;0];
bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(3)*AOAScale;2*pi];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;

zetaF = deg2rad(97);

alpha0 = 0;

% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF];
hf = 25000;
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF; hf; 0];	
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF; hf];	
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF; hf; 0];	

bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF];
%  bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF; 0];
bounds.upper.events = bounds.lower.events;

% bounds.lower.path = [49999];	
% bounds.upper.path = [50001];

%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

MoonLander.bounds = bounds;

%------------------------------------------------------
% Select the number of nodes for the spectral algorithm
%------------------------------------------------------

% algorithm.nodes = [60];  % somewhat arbitrary number; theoretically, the 
  algorithm.nodes = [90];                        % larger the number of nodes, the more accurate 
                         % the solution (but, practically, this is not
                         % always true!)
%70 for 55kPa, 90 for others


   
t0			= 0;
tfGuess 	= 90;			
% slightly educated guess of final time (for the scaled problem!)

guess.states(1,:)	= [h0, 25500]; %24.5 for 45kpa, 23 for 55kpa
guess.states(2,:)	= [v0, 1500];
guess.states(3,:)	= [m0, mF];
guess.states(4,:)	= [gamma0,0];
guess.states(5,:)	= [0, deg2rad(-2)];
guess.states(6,:)	= [1.63, zetaF];
% guess.states(5,:)	= [0, 0];
guess.controls		= [0.0, 0.0];
guess.time			= [t0, tfGuess];


%% Start Iterative Plot
figure(10)
plot(linspace(guess.time(1),guess.time(2),algorithm.nodes),linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes),'k')
iterative_V(end+1,:) = linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes);
iterative_t(end+1,:) = linspace(guess.time(1),guess.time(2),algorithm.nodes);

filename = 'testnew51.gif';
frame = getframe(10);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',0);
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
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [rad2deg(y(:,4).')/100 rad2deg(gamma)/100],'color','k','linestyle','-');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [y(:,2).'/1000 v/1000],'color','k','linestyle','--');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [y(:,1).'/10000 V/10000],'color','k','linestyle',':');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha)/10],'color','k','linestyle','-.');
legend('Trajectory Angle (degrees/100)','Velocity (m/s / 1000)','Altitude (km /10)','AoA (degrees/10)')
xlabel('Time (s)')
xlim([0,primal.nodes(end)+t_prepitch(end)]);

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

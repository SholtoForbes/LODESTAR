
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
clc
%==============================================================
% global zeta
global phi % latitude
global xi % longitude
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
MoonLander.cost 		= 'FirstStageReturnCost';
MoonLander.dynamics	    = 'FirstStageReturnDynamics';
MoonLander.events		= 'FirstStageReturnEvents';	
%MoonLander.path		= 'FirstStageReturnPath';
global scattered

addpath TrajOpt-master

Aero = dlmread('GlideBackAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

% M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
% M_interp = unique(sort(Aero(:,1)));
% 
% AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
% AoA_interp = unique(sort(Aero(:,2)));
% 
% [grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
% grid.Lift = scattered.Lift(grid.M,grid.AoA);
% grid.Drag = scattered.Drag(grid.M,grid.AoA);
% scattered.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
% scattered.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');

global SPARTANscale

SPARTANscale = 1;


% TARGET ==================================================================
%target final altitude and trajectory angle
gammaf = 0; % set final desired flight angle, JC change from 0
% =========================================================================

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
global AOAScale
AOAScale = 1;

%Define initial conditions
h0 = 25003;  
v0 = 1492.8;   
gamma0 = 0;    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

vf = 67;

hLow = 1300;   %
hUpp = 50000;  

vLow = 0; 
vUpp = 2500;  

gammaLow = deg2rad(-89.9);
gammaUpp = deg2rad(89.9);
% gammaLow = deg2rad(-45);
% gammaUpp = deg2rad(45);

zeta0 = 1.6;

% This sets the control limits, this is second derivative of AoA
% uLow = [-.001]*AOAScale; % Can do either AoA or thrust
% uUpp = [.001]*AOAScale; 
uLow = [-.01]*AOAScale; % Can do either AoA or thrust
uUpp = [.01]*AOAScale; 
% uLow = [-.1]*AOAScale; % Can do either AoA or thrust
% uUpp = [.1]*AOAScale;
%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 800;     % large upper bound; do not choose Inf

bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 

% These define the search space of the solution, including maximum AoA limits
bounds.lower.states = [hLow; vLow; gammaLow;-deg2rad(5)*AOAScale;0;-0.1];
bounds.upper.states = [ hUpp;  vUpp; gammaUpp;deg2rad(10)*AOAScale;2*pi; 0.1];
% bounds.lower.states = [hLow; vLow; gammaLow;-deg2rad(5)*AOAScale;0;-1];
% bounds.upper.states = [ hUpp;  vUpp; gammaUpp;deg2rad(10)*AOAScale;2*pi; 1];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;

% alpha0 = 0; %Set initial angle of attack to 0

% ========================================================================
 bounds.lower.events = [h0; v0; gamma0; zeta0; vf ; gammaf];
bounds.upper.events = bounds.lower.events;


% To be used if a 50kPa end state is desired
% bounds.lower.path = [49999];	
% bounds.upper.path = [50001];

%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

MoonLander.bounds = bounds;

%------------------------------------------------------
% Select the number of nodes for the spectral algorithm
%------------------------------------------------------

  % somewhat arbitrary number; theoretically, the 
                       % larger the number of nodes, the more accurate 
                         % the solution (but, practically, this is not
                         % always true!)
  algorithm.nodes = [100]; 
  % Change this by a few nodes to potentially change the solution slightly
  
  
  
% Initial Guess ===========================================================   
t0			= 0;
tfGuess 	= 300;	% JC changed from 90		
% slightly educated guess of final time (for the scaled problem!)

guess.states(1,:)	= [h0, h0]; %24.5 for 45kpa, 23 for 55kpa
guess.states(2,:)	= [v0, vf];
guess.states(3,:)	= [gamma0,0];
guess.states(4,:)	= [deg2rad(0), deg2rad(0)];
guess.states(5,:)	= [1.6, 1.597];

guess.states(6,:)	= [0, 0];

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
gamma = primal.states(3,:);
alpha = primal.states(4,:);
zeta = primal.states(5,:);
figure;
hold on
plot([primal.nodes], [ rad2deg(gamma)/100],'color','k','linestyle','-');
plot([primal.nodes], [ v/1000],'color','k','linestyle','--');
plot([primal.nodes], [ V/10000],'color','k','linestyle',':');
plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.');
legend('Trajectory Angle (degrees/100)','Velocity (m/s / 1000)','Altitude (km /10)','AoA (degrees/10)')
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
% [t_postpitch_f, postpitch_f] = ode45(@(t,postpitcqh_f) rocketDynamics(postpitch_f,ControlFunction(t,primal.nodes,dalphadt),phase,scattered), tspan, postpitch0_f);


% Forward Integrator
 phase = 'postpitch';
tspan = primal.nodes; 
postpitch0_f = [V(1) v(1) gamma(1) alpha(1) zeta(1) phi(1)];
[t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t,primal.nodes,alpha),phase,scattered), tspan, postpitch0_f);






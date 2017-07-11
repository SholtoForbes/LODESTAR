function [states_end] = FirstStageProblem(hf,gammaf,phif,zetaf,const)

global Throttle
Throttle = 0.85; % throttle the Merlin engine down by a constant value, to enable easier pitchover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;		
% clc
%==============================================================
% global zeta
% global phi
global q_forward

global iterative_V
iterative_V = [];
global iterative_t
iterative_t = [];
global iteration
iteration = 1;
global iterative_V_f
iterative_V_f = [];

global mach
%-----------------------------------
% Define the problem function files:
%-----------------------------------
MoonLander.cost 		= 'FirstStageCost';
MoonLander.dynamics	    = 'FirstStageDynamics';
MoonLander.events		= 'FirstStageEvents';	
% MoonLander.path		= 'LanderPath';
global scattered

addpath TrajOpt-master

Aero = dlmread('FirstStageAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));
% scattered.Moment = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = scattered.Lift(grid.M,grid.AoA);
grid.Drag = scattered.Drag(grid.M,grid.AoA);
% grid.Moment = scattered.Moment(grid.M,grid.AoA);
scattered.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
scattered.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');
% scattered.MomentGridded = griddedInterpolant(grid.M,grid.AoA,grid.Moment,'spline','linear');

global SPARTANscale

SPARTANscale = 1;


% TARGET ==================================================================
%target final altitude and trajectory angle
% global hf
% hf = 24000; % set final desired altitude
% hf = 24419; % const 1

% gammaf = 0; % set final desired flight angle

% gammaf =0.0408; % const 1
% =========================================================================



mRocket =21816 % total mass of scaled Falcon, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx
mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine) + mEngine; % structural mass fraction calculated without engine
mSpartan = 9819.11;


mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% h0_prepitch = 0;  %Rocket starts on the ground
% v0_prepitch = 0;  %Rocket starts stationary
% m0_prepitch = mTotal;  %Rocket starts full of fuel
% gamma0_prepitch = deg2rad(90);
% 
% phase = 'prepitch';
% tspan = [0 25]; % time to fly before pitchover (ie. straight up)
% 
% y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0];
% 
% % this performs a forward simulation before pitchover. The end results of
% % this are used as initial conditions for the optimiser. 
% [t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,scattered), tspan, y0);  


% FOR TESTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%this is a forward simulation for a set amount of time, at a constant AoA
% use this o test if your rocket can actually get close to where you need it
% phase = 'postpitch';
% Tratio = 1;
% tspan = [0 mFuel/156]; % flies for way too long
% postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) deg2rad(-1.4) 1.63];
% [t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,0,0,phase,scattered), tspan, postpitch0);
% 
% y
% postpitch
% postpitch(end,4)
% 

y = [90 30 0 0 0 0 0] % assume pitchover conditions, will need to back-simulate the necessary pitchover time
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
global AOAScale
AOAScale = 1;

%Define initial conditions
h0 = y(end,1);  
v0 = y(end,2);  
m0 = y(end,3);  
gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 
% gamma0 = deg2rad(89); 

mF = mEmpty+mSpartan;  %Assume that we use all of the fuel

hLow = 0;   %Cannot go through the earth
hUpp = 30000;  

vLow = 0; 
vUpp = 3000;  

mLow = mEmpty;
mUpp = mTotal;
% if const == 32
%     mUpp = 30000;
% else
% mUpp = 29000;
% end
% mUpp = 29000;
gammaLow = deg2rad(-.1);
% gammaUpp = deg2rad(89.9);
gammaUpp = gamma0;
% This sets the control limits, this is second derivative of AoA
% uLow = [-.001]*AOAScale; % Can do either AoA or thrust
% uUpp = [.001]*AOAScale; 
uLow = [-.0005]*AOAScale; % Can do either AoA or thrust
uUpp = [.0005]*AOAScale; 
%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 300;     % large upper bound; do not choose Inf

bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 

% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.

% bounds.lower.states = [hLow; vLow; mF-1;-0.1;-deg2rad(4)*AOAScale;0];
% bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(3)*AOAScale;2*pi];


% These define the search space of the solution, including maximum AoA limits
bounds.lower.states = [hLow; vLow; mF-1;gammaLow;-deg2rad(5)*AOAScale;0;-0.1; -0.25];
bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(2)*AOAScale;2*pi; 0.1; -0.15];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;

% zetaF = deg2rad(97); %Targets a specific heading angle. This is about right to get into a SSO by the end of third stage flight
% zetaF = 1.6837;
% phif = -0.2138; % these match second stage

alpha0 = 0; %Set initial angle of attack to 0


% Commented out bounds are for different end states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF];
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF; hf; 0];	
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; mF; zetaF; hf];	
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF; hf; 0];	
% bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF];
%   bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF; 0.0; hf];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  bounds.lower.events = [h0; v0; m0; gamma0; alpha0; zetaF; gammaf];
% bounds.lower.events = [h0; v0; gamma0; alpha0; zetaF; gammaf; 1500; mEmpty+mSpartan; hf; phif];
if const == 32
    vf = 1596;
else
vf = 1520;
end
Atmosphere = dlmread('atmosphere.txt');

bounds.lower.events = [h0; v0; gamma0; alpha0; zetaf; gammaf; vf; mEmpty+mSpartan; hf; phif];
% bounds.lower.events = [h0; v0; gamma0; alpha0; zetaf; gammaf; vf; hf; phif];

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
if const == 3
%     algorithm.nodes = [82]; 
algorithm.nodes = [90]; 
else
  algorithm.nodes = [82]; 
end
%  algorithm.nodes = [60]; 
  % Change this by a few nodes to potentially change the solution slightly
  
  
  
% Initial Guess ===========================================================   
t0			= 0;
% tfGuess 	= 90;			
tfGuess 	= 80;
% slightly educated guess of final time (for the scaled problem!)
% guess.states(1,:)	= [h0, hf+100]; %24.5 for 45kpa, 23 for 55kpa

% if const == 3 || const == 12
    guess.states(1,:)	= [h0, hf-100];
% else
%     guess.states(1,:)	= [h0, hf-500];
% end
% guess.states(1,:)	= [hf,hf];
% guess.states(1,:)	= [h0, interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/vf^2)+100]; %24.5 for 45kpa, 23 for 55kpa
guess.states(2,:)	= [v0, 1520];
% guess.states(2,:)	= [v0, 1620];
% guess.states(2,:)	= [v0, 1550];
% guess.states(3,:)	= [28000, mF];


% if const == 3
    guess.states(3,:)	= [28000, mF];
%     else
%     guess.states(3,:)	= [mTotal, mF];
% end

% guess.states(3,:)	= [mTotal, mF];
guess.states(4,:)	= [gamma0,gammaf];
% guess.states(5,:)	= [deg2rad(-1), deg2rad(-2)];
guess.states(5,:)	= [deg2rad(0), deg2rad(-2)];
guess.states(6,:)	= [1.63, zetaf];

guess.states(7,:)	= [0, 0];

guess.states(8,:)	= [-0.22, phif];

guess.controls		= [0.0, 0.0];
guess.time			= [t0, tfGuess];


%% Start Iterative Plot
figure(1010)
plot(linspace(guess.time(1),guess.time(2),algorithm.nodes),linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes),'k')
iterative_V(end+1,:) = linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes);
iterative_t(end+1,:) = linspace(guess.time(1),guess.time(2),algorithm.nodes);

filename = 'testnew51.gif';
frame = getframe(1010);
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
phi = primal.states(8,:);

%=============================================================================

% FOR TESTINGm, see where it gets 
% dalphadt = [diff(alpha)./diff(primal.nodes) 0];
% 
% phase = 'postpitch';
% Tratio = 1;
% tspan = primal.nodes; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) alpha(1)];
% [t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamics(postpitch_f,ControlFunction(t,primal.nodes,dalphadt),phase,scattered), tspan, postpitch0_f);


%% Forward Integrator
 phase = 'postpitch';
tspan = primal.nodes; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [y(end,1) y(end,2) m(1) deg2rad(89.9) phi(1) zeta(1)];

[t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t,primal.nodes,zeta),ControlFunction(t,primal.nodes,alpha),phase,scattered,Throttle), tspan, postpitch0_f);

figure(103)
hold on
plot(postpitch_f(:,1));
plot(V);
% y
% postpitch_f
% postpitch_f(end,4)


% fuel_left = mFuel - (m(1) - m(end)) % for fixed mass


% Iterative Prepitch Determination ========================================
%This back determines the mass and launch altitude necessary to get to
%100m, 30m/s at the PS method determined fuel mass

% ntoe that launch altitude does vary, but it should only be slightly
controls = fminunc(@(controls) prepitch(controls,m(1),scattered,Throttle),[10,6]);

h_launch = controls(1)
t_prepitch = controls(2)
Isp = 275;
T = 422581;
dm = -T./Isp./9.81;
m0_prepitch = m(1) - dm*t_prepitch;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = h_launch;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
% m0_prepitch = m0_prepitch;  %Rocket starts full of fuel
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan = [0 t_prepitch]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0, 0];

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,scattered,Throttle), tspan, y0);  


figure(101);
hold on
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [rad2deg(y(:,4).')/100 rad2deg(gamma)/100],'color','k','linestyle','-');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [y(:,2).'/1000 v/1000],'color','k','linestyle','--');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [y(:,1).'/10000 V/10000],'color','k','linestyle',':');
plot([t_prepitch.' primal.nodes+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha)/10],'color','k','linestyle','-.');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')

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

figure(102)

plot(t,dLHdu,t,mu_1,t,mu_2,t,mu_3,t,mu_4,t,mu_5,t,mu_u);
legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
title('Validation')


states_end = [t_prepitch.' primal.nodes+t_prepitch(end) ; y.' primal.states];

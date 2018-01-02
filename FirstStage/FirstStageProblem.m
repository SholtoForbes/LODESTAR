function [states_end] = FirstStageProblem(hf,gammaf,phif,zetaf,const)

global Throttle
Throttle = 0.85; % throttle the Merlin engine down by a constant value, to enable easier pitchover

% Aerodynamics File Path
Aero = dlmread('FirstStageAeroCoeffs.txt');


%% Vehicle 
global Vehicle 

mRocket =21816 % total mass of scaled Falcon, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx
mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine) + mEngine; % structural mass fraction calculated without engine
mSpartan = 9819.11;

% Thrust and Isp are modified with altitude through the formula:
% SL + (101325-P_atm)*Mod

Vehicle.T.SL = 555900; % Thrust from Falcon 1 users guide. 
Vehicle.T.Mod = 0.5518; % exit area calculated in SCALING.docx

Vehicle.Isp.SL = 275; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
Vehicle.Isp.Mod = 2.9410e-04;

global SPARTANscale
SPARTANscale = 1;
Vehicle.Area = 62.77*SPARTANscale^(2/3); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;		
% clc
%==============================================================
% global zeta
% global phi

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
FirstStage.cost 		= 'FirstStageCost';
FirstStage.dynamics	    = 'FirstStageDynamics';
FirstStage.events		= 'FirstStageEvents';	

global interp

%% Import Atmosphere
global Atmosphere
Atmosphere = dlmread('atmosphere.txt');

%% Calculate Aerodynamic Splines

interp.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
interp.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = interp.Lift(grid.M,grid.AoA);
grid.Drag = interp.Drag(grid.M,grid.AoA);
interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');

mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)

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

%% Assign Pitchover Conditions
y = [90 30 0 0 0 0 0] % assume pitchover conditions, launch conditions are back-simulated 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Define initial conditions
h0 = y(end,1);  
v0 = y(end,2);  
m0 = y(end,3);  
gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

mF = mEmpty+mSpartan;  %Assume that we use all of the fuel

hLow = 0;   %Cannot go through the earth
hUpp = 30000;  

vLow = 0; 
vUpp = 3000;  

mLow = mEmpty;
mUpp = mTotal;

phiL = -0.5;
phiU = -0.2;

if phif > phiU || phif < phiL
    disp('Given latitude outside of bounds')
end


gammaLow = deg2rad(-.1);
gammaUpp = gamma0;
% This sets the control limits, this is second derivative of AoA
uLow = [-.0005]; % Can do either AoA or thrust
uUpp = [.0005];

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax 	    = 300;     % large upper bound; do not choose Inf

bounds.lower.time 	= [0 0];				
bounds.upper.time	= [0 tfMax];			 

% Note: DO NOT set mfmin to zero because the equations 
% of motion have a singularity at m = 0.

% These define the search space of the solution, including maximum AoA limits
bounds.lower.states = [hLow; vLow; mF-1;gammaLow;-deg2rad(5);0;-0.1;phiL ];
bounds.upper.states = [ hUpp;  vUpp; mUpp;gammaUpp;deg2rad(2);2*pi; 0.1; phiU];

bounds.lower.controls = uLow;
bounds.upper.controls = uUpp;
alpha0 = 0; %Set initial angle of attack to 0


if const == 32
    vf = 1596;
else
vf = 1520;
end
Atmosphere = dlmread('atmosphere.txt');

bounds.lower.events = [h0; v0; gamma0; alpha0; zetaf; gammaf; vf; mEmpty+mSpartan; hf; phif];

bounds.upper.events = bounds.lower.events;


% To be used if a 50kPa end state is desired
% bounds.lower.path = [49999];	
% bounds.upper.path = [50001];

%------------------------------------
% Tell DIDO the bounds on the problem
%------------------------------------

FirstStage.bounds = bounds;

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
elseif const == 13
    algorithm.nodes = [91]; 
  elseif const == 32
    algorithm.nodes = [92];   
else
  algorithm.nodes = [90]; 
end
%  algorithm.nodes = [60]; 
  % Change this by a few nodes to potentially change the solution slightly
  
  
  
% Initial Guess ===========================================================   
t0			= 0;		
tfGuess 	= 80;
% educated guess of final time (for the scaled problem!)


guess.states(1,:)	= [h0, hf-100];
guess.states(2,:)	= [v0, 1520];
guess.states(3,:)	= [28000, mF];
guess.states(4,:)	= [gamma0,gammaf];
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
[cost, primal, dual] = dido(FirstStage, algorithm);
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

[t_postpitch_f, postpitch_f] = ode45(@(t,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t,primal.nodes,zeta),ControlFunction(t,primal.nodes,alpha),phase,interp,Throttle,Vehicle,Atmosphere), tspan, postpitch0_f);

figure(103)
hold on
plot(postpitch_f(:,1));
plot(V);

% Iterative Prepitch Determination ========================================
%This back determines the mass and launch altitude necessary to get to
%100m, 30m/s at the PS method determined fuel mass

% ntoe that launch altitude does vary, but it should only be slightly
controls = fminunc(@(controls) prepitch(controls,m(1),interp,Throttle,Vehicle,Atmosphere),[10,6]);

h_launch = controls(1)
t_prepitch = controls(2)
Isp = Vehicle.Isp.SL;
T = Vehicle.T.SL;
dm = -T./Isp./9.81;
m0_prepitch = m(1) - dm*t_prepitch;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = h_launch;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan = [0 t_prepitch]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0, 0];

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,interp,Throttle,Vehicle,Atmosphere), tspan, y0);  


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


clear all; clc

hf = 25000
gammaf = 0
phif = -0.4
zetaf = deg2rad(70)

const = 1



auxdata.Throttle = 0.85; % throttle the Merlin engine down by a constant value, to enable easier pitchover

% Aerodynamics File Path
Aero = dlmread('FirstStageAeroCoeffs.txt');


%% Vehicle 


mRocket =21816 % total mass of scaled Falcon, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx
mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine) + mEngine; % structural mass fraction calculated without engine
mSpartan = 9819.11;

% Thrust and Isp are modified with altitude through the formula:
% SL + (101325-P_atm)*Mod

auxdata.Vehicle.T.SL = 555900; % Thrust from Falcon 1 users guide. 
auxdata.Vehicle.T.Mod = 0.5518; % exit area calculated in SIZING.docx

auxdata.Vehicle.Isp.SL = 275; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
auxdata.Vehicle.Isp.Mod = 2.9410e-04;

auxdata.Vehicle.Area = 62.77; 


%% Import Atmosphere

auxdata.Atmosphere = dlmread('atmosphere.txt');

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

auxdata.interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
auxdata.interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');

mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)


%% Assign Pitchover Conditions

%Define initial conditions at pitchover, these are assumed
h0 = 90;  
v0 = 30;    

gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

mF = mEmpty+mSpartan;  %Assume that we use all of the fuel


alpha0 = 0; %Set initial angle of attack to 0

if const == 32
    vf = 1596;
else
vf = 1520;
end



hLow = 1;   %Cannot go through the earth
hUpp = 30000;  

vLow = 0; 
vUpp = 3000;  

mLow = mEmpty;
mUpp = mTotal;

phiL = -0.5;
phiU = -0.2;

zetaMin = 0;
zetaMax = 2*pi;

alphaMin = -deg2rad(5);
alphaMax = deg2rad(2);

dalphadt2Min = -0.1;
dalphadt2Max = 0.1;

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
		 

%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%


bounds.phase.initialtime.lower = 0;
bounds.phase.initialtime.upper = 0;

bounds.phase.finaltime.lower = 0;
bounds.phase.finaltime.upper = tfMax;

bounds.phase.initialstate.lower = [h0, v0,  mF-1, gamma0, alpha0,  zetaMin, dalphadt2Min, phiL ];
bounds.phase.initialstate.upper = [h0, v0, mUpp, gamma0, alpha0, zetaMax, dalphadt2Max, phiU];

bounds.phase.state.lower = [hLow, vLow, mF-1, gammaLow, alphaMin, zetaMin, dalphadt2Min, phiL ];
bounds.phase.state.upper = [ hUpp,  vUpp, mUpp, gammaUpp, alphaMax, zetaMax, dalphadt2Max, phiU];

bounds.phase.finalstate.lower = [hf, vf, mF, gammaf, alphaMin,  zetaf, dalphadt2Min,  phif];
bounds.phase.finalstate.upper = [hf, vf, mF, gammaf, alphaMax,  zetaf, dalphadt2Max,  phif];

bounds.phase.control.lower = uLow;
bounds.phase.control.upper = uUpp;

bounds.phase.path.lower = 0;
bounds.phase.path.upper = 50000;


guess.phase.time =  [0; tfMax];

guess.phase.state(:,1) = [h0; h0];
guess.phase.state(:,2) = [v0; vf];
guess.phase.state(:,3) = [mUpp; mF];
guess.phase.state(:,4) = [gamma0; gammaf];
guess.phase.state(:,5) = [alpha0; 0];
guess.phase.state(:,6) = [zetaf; zetaf];
guess.phase.state(:,7) = [0; 0];
guess.phase.state(:,8) = [phif; phif];

guess.phase.control = [0; 0];


%-------------------------------------------------------------------------%
%----------- Assemble All Information into Setup Structure ---------------%
%-------------------------------------------------------------------------%

setup.mesh.method       = 'hp-LiuRao-Legendre';
setup.mesh.maxiterations = 4;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 50;

setup.nlp.ipoptoptions.maxiterations = 300;


setup.name = 'Goddard-Rocket-Problem';
setup.functions.continuous = @FirstStageContinuous;
setup.functions.endpoint = @FirstStageEndpoint;
setup.nlp.solver = 'ipopt';
setup.bounds = bounds;
setup.guess = guess;
setup.auxdata = auxdata;
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
setup.scales.method = 'automatic-bounds';

setup.mesh.tolerance = 1e-6;
setup.method = 'RPM-Differentiation';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);






%============================================================================

t = output.result.solution.phase.time.';

V = output.result.solution.phase.state(:,1).';
v = output.result.solution.phase.state(:,2).';
m = output.result.solution.phase.state(:,3).';
gamma = output.result.solution.phase.state(:,4).';
alpha = output.result.solution.phase.state(:,5).';
zeta = output.result.solution.phase.state(:,6).';
phi = output.result.solution.phase.state(:,8).';


interp = auxdata.interp;
Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;
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
tspan = t; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [h0 v0 m(1) deg2rad(89.9) phi(1) zeta(1)];

[t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t_f,t,zeta),ControlFunction(t_f,t,alpha),phase,interp,Throttle,Vehicle,Atmosphere), tspan, postpitch0_f);

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
plot([t_prepitch.' t+t_prepitch(end)], [rad2deg(y(:,4).')/100 rad2deg(gamma)/100],'color','k','linestyle','-');
plot([t_prepitch.' t+t_prepitch(end)], [y(:,2).'/1000 v/1000],'color','k','linestyle','--');
plot([t_prepitch.' t+t_prepitch(end)], [y(:,1).'/10000 V/10000],'color','k','linestyle',':');
plot([t_prepitch.' t+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha)/10],'color','k','linestyle','-.');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')

legend('Trajectory Angle (degrees/100)','Velocity (m/s / 1000)','Altitude (km /10)','AoA (degrees/10)')
xlabel('Time (s)')
xlim([0,t(end)+t_prepitch(end)]);

% mu_1 = dual.states(1,:);
% mu_2 = dual.states(2,:);
% mu_3 = dual.states(3,:);
% mu_4 = dual.states(4,:);
% mu_5 = dual.states(5,:);
% 
% mu_u = dual.controls; % NOTE: This deviates from 0, as the controls are set as a buffer. Do not set a parameter directly tied to the vehicle model as the control.
% 
% %GRADIENT NORMALITY CONDITION
% 
% % Lagrangian of the Hamiltonian 
% dLHdu = dual.dynamics(3,:) + mu_u; % 

% figure(102)
% 
% plot(t,dLHdu,t,mu_1,t,mu_2,t,mu_3,t,mu_4,t,mu_5,t,mu_u);
% legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
% title('Validation')


% states_end = [t_prepitch.' t+t_prepitch(end) ; y.' primal.states];





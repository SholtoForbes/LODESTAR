function FirstStageOutput = FirstStageMain(hf,gammaf,phif,zetaf,mode,Timestamp)
% clear all; clc

% hf = 25000
% gammaf = 0
% phif = -0.4
% zetaf = deg2rad(70)
% 
% mode = 1

zetaf

% auxdata.Throttle = 0.75; % throttle the Merlin engine down by a modeant value, to enable easier pitchover
 auxdata.Throttle = .7
% Aerodynamics File Path
Aero = dlmread('FirstStageAeroCoeffs.txt');


%% Vehicle 


% mRocket =21816; % total mass of scaled Falcon at 9.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx

mRocket =19569; % total mass of scaled Falcon at 8.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx

mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine); % structural mass fraction calculated without engine
% mFuel = 0.85*mRocket; % structural mass fraction calculated without engine
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

if mode == 32
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

zetaMin = -2*pi;
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

bounds.phase.finalstate.lower = [hf, vLow, mF, gammaf, alphaMin,  zetaf, dalphadt2Min,  phif];
bounds.phase.finalstate.upper = [hf, vUpp, mF, gammaf, alphaMax,  zetaf, dalphadt2Max,  phif];

bounds.phase.finalstate.lower = [hLow, vLow, mF-1, gammaLow, alphaMin, zetaMin, dalphadt2Min, phiL ];
bounds.phase.finalstate.upper = [ hUpp,  vUpp, mUpp, gammaUpp, alphaMax, zetaMax, dalphadt2Max, phiU];

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
setup.mesh.maxiterations = 2;
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


figure(111);
hold on
title('First Stage Trajectory');
    fig = gcf;
set(fig,'Position',[200 0 850 600])
subplot(4,2,1)
hold on
title('Trajectory Angle (deg)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [rad2deg(y(:,4).') rad2deg(gamma)],'color','k');
subplot(4,2,2)
hold on
title('Velocity (m/s)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [y(:,2).' v],'color','k');
subplot(4,2,3)
hold on
title('Altitude (km)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [y(:,1).'/1000 V/1000],'color','k');
subplot(4,2,4)
hold on
title('Angle of Attack (deg)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha)],'color','k');
subplot(4,2,5)
hold on
title('Mass (kg)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [y(:,3).' m],'color','k');
subplot(4,2,6)
hold on
title('Heading Angle (deg)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [rad2deg(y(:,6).') rad2deg(zeta)],'color','k');
subplot(4,2,7)
hold on
title('Latitude (deg)');
xlim([0 t(end)+t_prepitch(end)]);
plot([t_prepitch.' t+t_prepitch(end)], [rad2deg(phi(1)+y(:,8).') rad2deg(phi)],'color','k');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')
xlabel('Time (s)')
xlim([0,t(end)+t_prepitch(end)]);

saveas(figure(111),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'FirstStage.fig']);

FirstStageOutput = output;


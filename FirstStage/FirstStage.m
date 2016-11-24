

clc; clear;
clear all;
addpath TrajOpt-master

Aero = dlmread('FirstStageAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));
   
mRocket = 27000; %(kg)  %Total lift-off mass
mFuel = 0.8*mRocket;  %(kg)  %mass of the fuel
mSpartan = 8755.1;
mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)
global Tmax
Tmax = 460000; %(nearly a merlin 1C) (This whole thing is nearly a Falcon 1 first stage)

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

% % FOR TESTING
% phase = 'postpitch';
% Tratio = .94;
% tspan = [0 (y(end,3)-(mEmpty+mSpartan))/(Tratio*60*Tmax/200000)];
% postpitch0 = [y(end,1), y(end,2), y(end,3), deg2rad(89)];
% [t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,Tratio*Tmax,phase), tspan, postpitch0);
% 
% y
% postpitch
% postpitch(end,4)
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

h0 = y(end,1);  %Rocket starts on the ground
v0 = y(end,2);  %Rocket starts stationary
m0 = y(end,3);  %Rocket starts full of fuel
gamma0 = deg2rad(88);    % pitchover 

vF = 1850;  
mF = mEmpty+mSpartan;  %Assume that we use all of the fuel
gammaF = deg2rad(1);
hF = 26550;
% hF = 30000;

hLow = 0;   %Cannot go through the earth
hUpp = 80000;  

vLow = 0; 
vUpp = inf;  

mLow = mEmpty;
mUpp = mTotal;

gammaLow = deg2rad(-1);
gammaUpp = deg2rad(90);

% uLow = [0];
% uUpp = [Tmax]; %Maximum thrust output

% uLow = [-1.5]; % Can do either AoA or thrust
% uUpp = [.5]; 

uLow = [-.1]; % Can do either AoA or thrust
uUpp = [.1]; 


P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;

P.bounds.finalTime.low = 0;
P.bounds.finalTime.upp = 60*60;

P.bounds.state.low = [hLow;vLow;mLow;gammaLow;-0.4];
P.bounds.state.upp = [hUpp;vUpp;mUpp;gammaUpp;0.4];

P.bounds.initialState.low = [h0;v0;m0;gamma0;0];
P.bounds.initialState.upp = [h0;v0;m0;gamma0;0];

% P.bounds.finalState.low = [hLow;vF;mF;gammaF];
% P.bounds.finalState.upp = [hUpp;vF;mF;gammaF];

P.bounds.finalState.low = [hF;vLow;mF;gammaF;0];
P.bounds.finalState.upp = [hF;vUpp;mF;gammaF;0];

P.bounds.control.low = uLow;
P.bounds.control.upp = uUpp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Initial Guess                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
hGuess = hF;   %(m) guess at the maximum height reached
P.guess.time = [0, 142];  %(s)
P.guess.state = [ [h0;v0;m0;gamma0;0],  [hGuess;vF;mF;gammaF;0] ];
P.guess.control = [ 0, 0 ];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Objective and Dynamic functions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Dynamics function:
phase = 'postpitch';
P.func.dynamics = @(t,x,u)( rocketDynamics(x,u,phase,scattered) );

% Objective function:
% P.func.bndObj = @(t0,x0,tF,xF)( -xF(1)/10000 );  %Maximize final height
P.func.bndObj = @(t0,x0,tF,xF)( -xF(2)/100 );
% P.func.bndObj = @(t0,x0,tF,xF)( 0 );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Options and Method selection                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


        P.options(1).method = 'rungeKutta';
        P.options(1).defaultAccuracy = 'medium';
        P.options(1).nlpOpt.MaxFunEvals = 2e5;
        P.options(1).nlpOpt.MaxIter = 1e5;


        
        P.options(2).method = 'rungeKutta';
        P.options(2).defaultAccuracy = 'high';
        P.options(2).nlpOpt.MaxFunEvals = 5e4;
        P.options(2).nlpOpt.MaxIter = 1e5;
        P.options(2).rungeKutta.nSegment = 40;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = trajOpt(P);

t = linspace(soln(end).grid.time(1),soln(end).grid.time(end),100);  % It interpolates the end result!
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);








% Forward Simulation ======================================================
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% %                        Pre-Pitchover Simulation                         %
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% % Note this doesnt work for AoA control
% 
% f_h0_prepitch = 0;  %Rocket starts on the ground
% f_v0_prepitch = 0;  %Rocket starts stationary
% f_m0_prepitch = mTotal;  %Rocket starts full of fuel
% f_gamma0_prepitch = deg2rad(90);
% 
% phase = 'prepitch';
% f_tspan = [0 15];
% f_y0 = [f_h0_prepitch, f_v0_prepitch, f_m0_prepitch, f_gamma0_prepitch, 0];
% % [f_t_prepitch, f_y_prepitch] = ode45(@(f_t,f_y) rocketDynamics(f_y,Tmax,phase), f_tspan, f_y0);
% [f_t_prepitch, f_y_prepitch] = ode45(@(f_t,f_y) rocketDynamics(f_y,0,phase), f_tspan, f_y0);
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% %                        Post-Pitchover Simulation                         %
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% f_h0 = f_y_prepitch(end,1);  %Rocket starts on the ground
% f_v0 = f_y_prepitch(end,2);  %Rocket starts stationary
% f_m0 = f_y_prepitch(end,3);  %Rocket starts full of fuel
% f_gamma0 = deg2rad(89);    % pitchover 
% 
% phase = 'postpitch';
% f_tspan = [0 t(end)];
% f_y0 = [f_h0, f_v0, f_m0, f_gamma0, 0];
% [f_t, f_y] = ode45(@(f_t,f_y) rocketDynamics(f_y,ControlFunction(f_t,t,u),phase), f_tspan, f_y0);




% Plotting

figure(120);
subplot(2,3,1);
hold on
plot(t,x(1,:)/1000)
% plot(f_t,f_y(:,1)/1000)
xlabel('time (s)')
ylabel('height (km)')
subplot(2,3,2);
plot(t,x(3,:))
xlabel('time (s)')
ylabel('mass (kg)')
subplot(2,3,3);
plot(t,x(2,:))
xlabel('time (s)')
ylabel('velocity (m/s)')
subplot(2,3,4);
plot(t,x(4,:))
xlabel('time (s)')
ylabel('trajectory angle (rad)')
subplot(2,3,5);
plot(t,u)
xlabel('time (s)')
ylabel('Control')
subplot(2,3,6);
plot(t,x(5,:))
xlabel('time (s)')
ylabel('Alpha (rad)')

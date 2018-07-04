


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
clc
%==============================================================
	

global scattered

addpath TrajOpt-master

Aero = dlmread('FirstStageAeroCoeffs.txt');
scattered.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
scattered.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

global SPARTANscale
% SPARTANscale = 0.75
SPARTANscale = 1

% mRocket = 27000; %(kg)  %Total lift-off mass
mRocket = 40000; %(kg)  %Total lift-off mass
mFuel = 0.8*mRocket;  %(kg)  %mass of the fuel
% mFuel = 0.7814*mRocket-630*1.1^(3/2);  %(kg)  %mass of the fuel   630 is for the merlin , change this if scaling
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

% FOR TESTINGm, see where it gets after 100s no aoa
phase = 'postpitch';
Tratio = 1;
% tspan = [0 mFuel/156];
tspan = [0 200];
postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) deg2rad(0.05)];
[t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,0,phase,scattered), tspan, postpitch0);

y
postpitch
postpitch(end,4)



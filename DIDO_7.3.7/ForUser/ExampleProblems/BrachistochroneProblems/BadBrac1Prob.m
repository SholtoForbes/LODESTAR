%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem (script) file for the Bad Brachistochrone Problem Brac: 1
% Template for A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;				% always a good idea to begin with this!	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============================
% Problem variables:
%----------------------------
% states = (xBar, yBar, vBar)
% controls = thetaBar
%============================

%------------------------------------------
% Globalize the Problem Constants and Units
%------------------------------------------

global CONSTANTS UNITS SCALES

%==================================
% Load constants, units and scales.
%==================================

%----------------------------------------------------------------------
% In a complex problem, you may want to load such quantites from
% a separate file so that this "problem" file is clean and modular.
%----------------------------------------------------------------------

%-------------------------------------------
CONSTANTS.g = 9.81;      % Earth g in m/s^2
%-------------------------------------------
%---------------
% Bad Brac data:
%---------------
xFinal = 1000;     % in meters
yFinal = 1;        % in meters
%-----------------------------
%--------------------------------------------------------------------------
UNITS.x = 100;      % x unit is based on xFinal
UNITS.y = 20;       % y unit is based on DIDO runs for good values of xFinal
UNITS.v = 10;       % ditto
UNITS.theta = 1;
UNITS.t = 10;
% %--------------------------------------------------------------------------
% % Good data for testing:
% xFinal = 10;        % See Brac 1 Problem File
% yFinal = 10;
% UNITS.x = 1;      
% UNITS.y = 1;       
% UNITS.v = sqrt(9.81);       
% UNITS.theta = 1;
% UNITS.t = sqrt(1/9.81);
% %-----------------------------------------
SCALES.xdot = UNITS.t*UNITS.v/UNITS.x;
SCALES.ydot = UNITS.t*UNITS.v/UNITS.y;
SCALES.vdot = UNITS.t/UNITS.v;
%------------------------------------------

%---------------------------------------
% bounds the state and control variables
%---------------------------------------

xLBar = 0;  xUBar = 10*UNITS.x;
yLBar = 0;  yUBar = 10*UNITS.y;
vLBar = 0;  vUBar = 100*UNITS.v;

bounds.lower.states = [xLBar; yLBar; vLBar];
bounds.upper.states = [xUBar; yUBar; vUBar];

bounds.lower.controls = [0];
bounds.upper.controls = [pi/UNITS.theta];


%------------------
% bound the horizon
%------------------
t0	    = 0;
tUBar 	= 200/UNITS.t;                           % swag

bounds.lower.time 	= [t0 t0];				
bounds.upper.time	= [t0 tUBar];			    % Fixed time at t0 and a possibly free time at tf


%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------
% See events file for definition of events function

x0Bar = 0; y0Bar = 0; v0Bar = 0;

bounds.lower.events = [x0Bar, y0Bar, v0Bar, xFinal/UNITS.x, yFinal/UNITS.y].';
bounds.upper.events = bounds.lower.events;      % equality event function bounds


%===========================================
% Define the problem using DIDO expressions:
%===========================================
Brac_1.cost 		= 'Brac1Cost';
Brac_1.dynamics	    = 'BadBrac1Dynamics';
Brac_1.events		= 'BadBrac1Events';		
%Path file not required for this problem formulation;

Brac_1.bounds       = bounds;
%====================================================

algorithm.nodes		= [24];					    % represents some measure of desired solution accuracy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call dido
tStart= cputime;    % start CPU clock 
[cost, primal, dual] = dido(Brac_1, algorithm);
runTime = cputime-tStart
% Ta da!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============================
% Debugging plots:
%=============================
xBar = primal.states(1,:);
yBar = primal.states(2,:);
vBar = primal.states(3,:);
thetaBar = primal.controls;
tBar = primal.nodes;

%---------------------------------------
figure;
plot(tBar, [xBar; yBar; vBar; thetaBar])
legend('Scaled States and Controls')
%----------------------------------------
figure;
plot(tBar, dual.dynamics);
legend('Scaled Costates');
%---------------------------------------
figure;
plot(tBar, dual.Hamiltonian);
legend('Scaled Hamiltonian');
%---------------------------------------

%=======================
% Desired plots:
%=======================

x = UNITS.x*primal.states(1,:);
y = UNITS.y*primal.states(2,:);
v = UNITS.v*primal.states(3,:);
t = UNITS.t*primal.nodes;

lamx = (UNITS.t/UNITS.x)*dual.dynamics(1,:);
lamy = (UNITS.t/UNITS.y)*dual.dynamics(2,:);
lamv = (UNITS.t/UNITS.v)*dual.dynamics(3,:);

OptimalCost = cost*UNITS.t


%----------------------------
% linear theoretical values:
%---------------------------
slopeAngle = atan(yFinal/xFinal);
LinearCost = sqrt(2*yFinal/(CONSTANTS.g*sin(slopeAngle)^2))

yLinear = yFinal/xFinal*x;
%============================================================================
figure;
plot(x, -y, '-o', x, -yLinear);
legend('DIDO Solution', 'Linear Solution');
title('Bad Brachistochrone Problem')
xlabel('x (m)');
ylabel('y (m)');

%--------------------------------------------------------------------------
figure;
plot(t, x, '-o', t, y, '-x', t, v, '-.');
title('"Bad" Brachistochrone States: Brac 1')
xlabel('time');
ylabel('states');
legend('x', 'y', 'v');

%----------------------------------------------
figure;
plot(t, dual.Hamiltonian);
title('"Bad" Brachistochrone Hamiltonian Evolution');
legend('H');
xlabel('time');
ylabel('Hamiltonian Value');

%----------------------------------------------
figure;
plot(t, [lamx; lamy; lamv]);
title('"Bad" Brachistochrone Costates: Brac 1')
xlabel('time');
ylabel('costates');
legend('\lambda_x', '\lambda_y', '\lambda_v');
%==============================================

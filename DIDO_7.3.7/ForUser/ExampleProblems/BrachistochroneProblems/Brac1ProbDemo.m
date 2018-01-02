%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem (script) file for the Brachistochrone Problem Formulation, Brac: 1
% Template for A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;				% always a good idea to begin with this!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================
% Problem variables:
%-------------------
% states = (x, y, v)
% controls = theta
%===================

%---------------------------------------
% bounds the state and control variables
%---------------------------------------

xL = 0; xU = 20;
yL = 0; yU = 20;
vL = 0; vU = 20;

bounds.lower.states = [xL; yL; vL];
bounds.upper.states = [xU; yU; vU];

bounds.lower.controls = [0];
bounds.upper.controls = [pi];


%------------------
% bound the horizon
%------------------
t0	    = 0;
tfMax 	= 10;   % swag for max tf; DO NOT set to Inf even for time-free problems

bounds.lower.time 	= [t0; t0];				
bounds.upper.time	= [t0; tfMax];			    % Fixed time at t0 and a possibly free time at tf


%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------
% See events file for definition of events function

bounds.lower.events = [0; 0; 0; 10; 10];
bounds.upper.events = bounds.lower.events;      % equality event function bounds


%============================================
% Define the problem using DIDO expresssions:
%============================================
Brac_1.cost 		= 'Brac1Cost';
Brac_1.dynamics	    = 'Brac1Dynamics';
Brac_1.events		= 'Brac1Events';		
%Path file not required for this problem formulation;	

Brac_1.bounds       = bounds;
%====================================================

algorithm.nodes		= [16];					    % represents some measure of desired solution accuracy

% guess.states = [0, 10;
%                 0, 10;
%                 0, 10];
% guess.controls = [0, pi/2];
% guess.time = [0, 2];

%load Brac1Data
%algorithm.guess = primal; 


% Call dido
tStart= cputime;    % start CPU clock 
[cost, primal, dual] = dido(Brac_1, algorithm);
runTime = cputime-tStart
% Ta da!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save Brac1Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = primal.states(1,:);
y = primal.states(2,:);
v = primal.states(3,:);
t = primal.nodes;

%=============================================
figure;
plot(t, x, '-o', t, y, '-x', t, v, '-.');
title('Brachistochrone States: Brac 1')
xlabel('time');
ylabel('states');
legend('x', 'y', 'v');

%----------------------------------------------
figure;
plot(t, dual.Hamiltonian);
title('Brachistochrone Hamiltonian Evolution');
legend('H');
xlabel('time');
ylabel('Hamiltonian Value');

%----------------------------------------------
figure;
plot(t, dual.dynamics);
title('Brachistochrone Costates: Brac 1')
xlabel('time');
ylabel('costates');
legend('\lambda_x', '\lambda_y', '\lambda_v');
%==============================================

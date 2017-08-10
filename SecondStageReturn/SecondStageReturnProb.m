%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SPARTAN_SCALE
SPARTAN_SCALE = 1 % volumetric scale

global scattered
% Counts iterations of DIDO
global iterative_V
iterative_V = [];
global iterative_t
iterative_t = [];
global iteration
iteration = 1;
global iterative_V_f
iterative_V_f = [];

%%
% Copy the current setting to archive
% This saves the entire problem file every time the program is run. 

Timestamp = datestr(now,30)

%% Inputs ============================================
%Take inputs of communicator matrices, these should be .txt files 
communicator = importdata('communicator.txt');
communicator_trim = importdata('communicator_trim.txt');


scattered.flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
scattered.flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
scattered.flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,6));


[MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';

scattered.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
scattered.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
scattered.pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');

% Produce Atmosphere Data
global Atmosphere
Atmosphere = dlmread('atmosphere.txt');



%=============================================== 

% V0 = 20000.;
% Vf = 50000.; %

%===================
% Problem variables:
% control factor: omega
% variable trajectory
%===================
 
%========================================================

%---------------------------------------
% bound and scale the state and control variables
%---------------------------------------

% Bounding the state and control variables creates a 'search space'. The
% optimal solution is found within this search space.These bounds should be sufficiently wide to allow a wide
% search space, but tight enough that a solution can be found efficiently.  These bounds must be
% chosen carefully, with the physical model in mind. 

VL = 10;
VU = 50000; 

% vL = 1500;
vL = 10;
vU = 1500; 


gammaL = -1;
gammaU = 1; % 



global scale

scale.V = 1;
scale.v = 1;
scale.gamma = 1;
scale.gammadot = 1;

alphaL = 0;
alphaU = deg2rad(20);

zetaL = 4;
zetaU = 6;

phiL = -.5;
phiU = -0.0;

xiL = -0.5;
xiU = 0.5;

% etaL = -.5;
% etaU = .5;
% etaL = -1;
% etaU = 1;

alphadotL = -0.01;
alphadotU = 0.01;

% etadotL = -0.01;
% etadotU = 0.01;

% 
% bounds.lower.states = [VL ; vL; gammaL; alphaL; zetaL; phiL; xiL; etaL];
% bounds.upper.states = [VU ; vU; gammaU; alphaU; zetaU; phiU; xiU; etaU];


bounds.lower.states = [VL ; vL; gammaL; alphaL; zetaL; phiL; xiL];
bounds.upper.states = [VU ; vU; gammaU; alphaU; zetaU; phiU; xiU];

% control bounds
 bounds.lower.controls = [alphadotL];
bounds.upper.controls = [alphadotU]; 

%  bounds.lower.controls = [alphadotL;etadotL];
% bounds.upper.controls = [alphadotU;etadotU]; 


%------------------
% bound the horizon
%------------------
% time bounds, this is unscaled

t0	    = 0;

 tfMax 	= 10000;

bounds.lower.time 	= [t0; 100];	
bounds.upper.time	= [t0; tfMax];

%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------

v0 = 1000;
% vf = 1550;


%% Define Events
V0 = 35000;
gamma0 = 0.048;
zeta0 = 4.7;
phi0 = -0.05;
xi0 = 0;
% zetaf = 1.6915;
% zetaf = 4.7124;


% bounds.lower.events = [V0;v0; gamma0;zeta0;phi0;xi0;0;deg2rad(-45)];
bounds.lower.events = [V0;v0; gamma0;zeta0;phi0;xi0;0];

% bounds.upper.events = [V0;v0; gamma0;zeta0;phi0;xi0;20000;0];
% bounds.upper.events = [V0;v0; gamma0;zeta0;phi0;xi0;0; 10000];
bounds.upper.events = bounds.lower.events;      % equality event function bounds

    bounds.lower.path = 0;
bounds.upper.path = 50000;



%% 
%============================================
% Define the problem using DIDO expresssions:
%============================================
% Call the files whih DIDO uses
TwoStage2d.cost 		= 'SecondStageReturnCost';
TwoStage2d.dynamics	    = 'SecondStageReturnDynamics';
TwoStage2d.events		= 'SecondStageReturnEvents';	
TwoStage2d.path		= 'SecondStageReturnPath';
TwoStage2d.bounds       = bounds;


%% Node Definition ====================================================
% The number of nodes is important to the problem solution. 
% While any number should theoretically work, in practice the choice of
% node number can have a large effect on results.
 

algorithm.nodes		= [150]; 
global nodes
nodes = algorithm.nodes;


%%  Guess =================================================================


% guess.states(1,:) = [35000 ,35000 ]; % test for new interpolation
% guess.states(1,:) = [V0 ,V0 ];
guess.states(1,:) = [V0 ,1000 ];
guess.states(2,:) = [v0, 100];

guess.states(3,:) = [0.05,0.00];

guess.states(4,:) = [deg2rad(19),deg2rad(19)];

guess.states(5,:) = [4.7,4.761];

guess.states(6,:) = [-0.1,-.14];

guess.states(7,:) = [0,-0.1];

% guess.states(8,:) = [1,1];

% guess.states(8,:) = [0,0];

% Control guess.
guess.controls(1,:)    = [0,0]; 
% guess.controls(2,:)    = [0,0]; 

guess.time        = [t0 ,550];
% Tell DIDO the guess
%========================
algorithm.guess = guess;
%=====================================================================================


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





%% Call DIDO
% This is where all the work is done
% =====================================================================
tStart= cputime;    % start CPU clock 
[cost, primal, dual] = dido(TwoStage2d, algorithm);

%%

% runTime = (cputime-tStart)
% runTime/60

EndTime = datestr(now,30)


V = primal.states(1,:)*scale.V;

v = primal.states(2,:)*scale.v;

t = primal.nodes;

gamma = primal.states(3,:);

Alpha = primal.states(4,:);
% Alpha = deg2rad(7)*ones(1,length(V));

zeta = primal.states(5,:);

phi = primal.states(6,:);

xi = primal.states(7,:);

% eta = rad2deg(primal.states(8,:));

alphadot = primal.controls(1,:)*scale.gamma;

% ===================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dfuel
dfuel


global M
global q
global Fd
global Fueldt
global lift

global Thrust
global flapdeflection

global a
global eq


% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(t(i+1) - t(i))*cos(gamma(i)) + H(i);
end


figure(201)

subplot(5,5,[1,10])
hold on
plot(H, V)
plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')


subplot(5,5,11)
hold on
plot(t, v)
plot(t(algorithm.nodes(1)), v(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Velocity (m/s)')


subplot(5,5,12)
plot(t, M)
title('Mach no')

subplot(5,5,13)
plot(t, q)
title('Dynamic Pressure (pa)')

subplot(5,5,14)
hold on
plot(t, rad2deg(gamma))
plot(t(algorithm.nodes(1)), rad2deg(gamma(algorithm.nodes(1))), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory Angle (Deg)')



subplot(5,5,15)
plot(t, Fd)
title('Drag Force')


% subplot(5,5,20)
% plot(t, flapdeflection)
% title('Flap Deflection (deg)')

subplot(5,5,21)
plot(t, Alpha)
title('Angle of Attack (deg)')

subplot(5,5,22);
plot(t, dual.dynamics);
title('costates')
xlabel('time');
ylabel('costates');
legend('\lambda_1', '\lambda_2', '\lambda_3');

subplot(5,5,23)
Hamiltonian = dual.Hamiltonian(1,:);
plot(t,Hamiltonian);
title('Hamiltonian')



subplot(5,5,25)
hold on
plot(t, rad2deg(alphadot))
title('Omegadot Control (Deg/s2)')


figure(202)
sp1 = subplot(2,6,[1,6]);
ax1 = gca; % current axes
hold on
plot(H/1000, V/1000,'Color','k')

title('Trajectory')
xlabel('Earth Normal Distance Flown (km)')
ylabel('Vertical Position (km)')

for i = 1:floor(t(end)/30)
    [j,k] = min(abs(t-30*i));
    str = strcat(num2str(round(t(k))), 's');
    text(H(k)/1000,V(k)/1000,str,'VerticalAlignment','top', 'FontSize', 10);
    
    plot(H(k)/1000, V(k)/1000, '+', 'MarkerSize', 10, 'MarkerEdgeColor','k')
end

plot(H(end)/1000, V(end)/1000, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','k')


dim = [.65 .45 .2 .2];

thirdstageexample_H = [0+H(end) (H(end)-H(end - 1))+H(end) 20*(H(end)-H(end - 1))+H(end) 40*(H(end)-H(end - 1))+H(end) 60*(H(end)-H(end - 1))+H(end) 80*(H(end)-H(end - 1))+H(end)]/1000; %makes a small sample portion of an arbitrary third stage trajectory for example
thirdstageexample_V = [0+V(end) (V(end)-V(end - 1))+V(end) 20*((V(end)-V(end -1)))+V(end) 40*((V(end)-V(end -1)))+V(end) 60*((V(end)-V(end -1)))+V(end) 80*((V(end)-V(end -1)))+V(end)]/1000;
plot(thirdstageexample_H, thirdstageexample_V, 'LineStyle', '--','Color','k');

hold on
sp2 = subplot(2,6,[7,9]);
xlabel('time (s)')

hold on
ax2 = gca; % current axes
xlim([min(t) max(t)]);

line(t, rad2deg(gamma),'Parent',ax2,'Color','k', 'LineStyle','-')

line(t, M,'Parent',ax2,'Color','k', 'LineStyle','--')

line(t, v./(10^3),'Parent',ax2,'Color','k', 'LineStyle','-.')

line(t, q./(10^4),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% line(t, heating_rate./(10^5),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% line(t, Q./(10^7),'Parent',ax1,'Color','k', 'LineStyle','-', 'lineWidth', 2.0)

% legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)',  'Q (Mj x 10)')
h = legend(ax2,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
rect1 = [0.12, 0.35, .25, .25];
set(h, 'Position', rect1)


sp3 = subplot(2,6,[10,12]);
xlabel('time (s)')
ax3 = gca;
xlim([min(t) max(t)]);
line(t, Alpha,'Parent',ax3,'Color','k', 'LineStyle','-')
% line(t, eta,'Parent',ax3,'Color','k', 'LineStyle','--')
line(t, zeta,'Parent',ax3,'Color','k', 'LineStyle',':')
% line(t, flapdeflection,'Parent',ax3,'Color','k', 'LineStyle','--')

% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax3, 'AoA (degrees)','heading Angle (degrees)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

% saveas(figure(202),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'SecondStage.fig']);

% figure(203)
% 
% subplot(2,5,[1,5]);
% 
% line(t, dual.dynamics(1,:),'Color','k', 'LineStyle','-');
% line(t, dual.dynamics(2,:),'Color','k', 'LineStyle','--');
% line(t, dual.dynamics(3,:),'Color','k', 'LineStyle','-.');
% line(t, dual.dynamics(4,:),'Color','k', 'LineStyle',':');
% line(t, dual.dynamics(5,:),'Color','k', 'LineStyle','-','LineWidth',2);
% title('costates')
% xlabel('time');
% ylabel('Costates');
% % axis([0,t(end),-1,1])
% legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4','\lambda_5');
% 
% subplot(2,5,[6,10])
% Hamiltonian = dual.Hamiltonian(1,:);
% plot(t,Hamiltonian,'Color','k');
% axis([0,t(end),-1,1])
% title('Hamiltonian')
% 
% 
% % save results
% dlmwrite('primal.txt', [primal.states;primal.controls;primal.nodes;q;IspNet;Alpha;M;eq;flapdeflection;phi]);
% dlmwrite('payload.txt', ThirdStagePayloadMass);
% dlmwrite('dual.txt', [dual.dynamics;dual.Hamiltonian]);
% dlmwrite('ThirdStage.txt',[ThirdStageZeta;ThirdStagePhi;ThirdStageAlt;ThirdStagev;ThirdStaget;[ThirdStageAlpha 0];ThirdStagem;ThirdStagegamma;[ThirdStageq 0]]);
% 
% 
% 
% copyfile('primal.txt',sprintf('../ArchivedResults/%s/primal_%s.txt',Timestamp,Timestamp))
% copyfile('dual.txt',sprintf('../ArchivedResults/%s/dual_%s.txt',Timestamp,Timestamp))
% copyfile('payload.txt',sprintf('../ArchivedResults/%s/payload_%s.txt',Timestamp,Timestamp))
% copyfile('ThirdStage.txt',sprintf('../ArchivedResults/%s/ThirdStage_%s.txt',Timestamp,Timestamp))
% primal_old = primal;
% 
% ts = timeseries(Isp,t);
% Mean_Isp = mean(ts)
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % TESTING AND VALIDATION
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % If these are valid then solution is a KKT point
% 
% %COMPLEMENTARY CONDITIONS
% % These should be zero if the state or control is within set bounds
% % <=0 if at min bound, >=0 if at max bound
% 
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
% 
% figure(205)
% 
% plot(t,dLHdu,t,mu_1,t,mu_2,t,mu_3,t,mu_4,t,mu_5,t,mu_u);
% legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
% title('Validation')
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % FORWARD INTEGRATION
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % This simply tests that the system dynamics hold, as the
% % Pseudospectral method may not converge to a realistic
% % solution
% 
% 
% 
% gamma_F = cumtrapz(t,gammadot)+ gamma(1);
% 
% gammadot_F = cumtrapz(t,omegadot) + gammadot(1);
% 
% 
% v_F = cumtrapz(t,a);
% v_F = v_F + v(1);
% 
% V_F = cumtrapz(t,v_F.*sin(gamma_F));
% V_F = V_F + V(1);
% 
% mfuel_F = cumtrapz(t,-Fueldt);
% mfuel_F = mfuel_F + mfuel(1);
% 
% figure(206)
% 
% subplot(5,1,1)
% plot(t,gamma_F,t,gamma);
% title('Forward Simulation Comparison');
% % note this is just a trapezoidal rule check, may not be exactly accurate
% subplot(5,1,2)
% plot(t,v_F,t,v);
% subplot(5,1,3)
% plot(t,V_F,t,V);
% subplot(5,1,4)
% plot(t,mfuel_F,t,mfuel);
% subplot(5,1,4)
% plot(t,gammadot_F,t,gammadot);
% 
% % Compute difference with CADAC for constant dynamic pressure path
% 
% t_diff = t - [0 t(1:end-1)];
% if const == 3
%     CADAC_DATA = dlmread('TRAJ.ASC');
%     CADAC_Alpha = interp1(CADAC_DATA(:,2),CADAC_DATA(:,4),M(1:68)); % 1:68 gives mach numbers that align, may need to change this
%     CADAC_V = interp1(CADAC_DATA(:,2),CADAC_DATA(:,11),M(1:68));
%     MeanError_V = sum(abs(CADAC_V - V(1:68))./V(1:68).*t_diff(1:68))/t(end)
%     MeanError_Alpha = sum(abs(CADAC_Alpha - Alpha(1:68))./Alpha(1:68).*t_diff(1:68))/t(end)
% end
% 
% % if PayloadGrid(phi(end),zeta(end),V(end)+10,gamma(end),v(end)) - PayloadGrid(phi(end),zeta(end),V(end),gamma(end),v(end)) < 0
% %     disp('Check Third Stage Payload Matrix, May Have Found False Maxima')
% % end
% if PayloadGrid(V(end)+10,gamma(end),v(end)) - PayloadGrid(V(end),gamma(end),v(end)) < 0
%     disp('Check Third Stage Payload Matrix, Found Maxima')
% end
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % FORWARD SIMULATION
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % This is a full forward simulation, using the angle of attack and flap
% % deflection at each node.
% 
% % Note, because the nodes are spaced widely, small interpolation
% % differences result in the forward simulation being slightly different
% % than the actual. This is mostly a check to see if they are close. 

mstruct = 4910.5 - 132.8 + 179.41; % mass of everything but fuel from dawids work
eta = 0;
m = mstruct;
forward0 = [V(1),gamma(1),v(1),zeta(1),phi(1),xi(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y, nodes,scattered, Atmosphere,ControlInterp(t,Alpha,f_t),eta),t(1:end),forward0);


figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(t,V);
% 
% 
% %% plot engine interpolation visualiser
% T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), V); 
% T1 = scattered.tempgridded(M,Alpha).*T0;
% M1 = scattered.M1gridded(M, Alpha);
% 
% plotM = [min(M_englist):0.01:9.5];
% plotT = [min(T_englist):1:550];
% [gridM,gridT] =  ndgrid(plotM,plotT);
% interpeq = scattered.eqGridded(gridM,gridT);
% interpIsp = scattered.IspGridded(gridM,gridT);
% 
% figure(210)
% hold on
% contourf(gridM,gridT,interpeq);
% scatter(data(:,1),data(:,2),30,data(:,4),'filled');
% plot(M1,T1,'r');
% 
% error_Isp = scattered.IspGridded(data(:,1),data(:,2))-data(:,3);
% 
% figure(211)
% hold on
% contourf(gridM,gridT,interpIsp);
% scatter(data(:,1),data(:,2),30,data(:,3),'filled')
% plot(M1,T1,'r');
% 
% %%
% [gridM2,gridAoA2] =  ndgrid(plotM,plotT);

%%

% =========================================================================
% Troubleshooting Procedure
% =========================================================================

% 1: Check that you have posed your problem correctly ie. it is physically
% feasible and the bounds allow for a solution
% 2: Check for NaN values (check derivatives in Dynamics file while running)
% 3: Check guess, is it reasonable? Is it too close to the expected
% solution? Both can cause errors! Sometimes there is no real rhyme or
% reason to picking the correct guess, but a close bound to
% the expected solution has worked the most in my experience
% 4: Play with the no. of nodes, try both even and odd values
% 5: Play with scaling
% 6: Try all of the above in various combinations until it works!






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;		
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
% Copy the current setting to archive
% This saves the entire problem file every time the program is run. 
% Disable this for conservatin of hard drive space
Timestamp = datestr(now,30)
copyfile('SecondStageProb.m',sprintf('../ArchivedResults/SecondStageProb_%s.m',Timestamp))
copyfile('SecondStageCost.m',sprintf('../ArchivedResults/SecondStageCost_%s.m',Timestamp))

% =========================================================================
% SET RUN MODE
% =========================================================================
% Change const to set the target of the simulation. Much of the problem
% definition changes with const.

% const = 1x: No end constraint, used for optimal trajectory calculation
% const = 1: 50kPa limit, 12: 55 kPa limit, 13: 45 kPa limit, 14: 50kPa limit & 10% additional drag

% const = 3: Fuel mass is constrained at end point, used for constant
% dynamic pressure calculation (50kPa constrained)
% const = 31: simple model for guess calc 

global const
const = 1

% Inputs ============================================
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

%produce scattered interpolants for thrust and fuel usage
% enginedata = dlmread('engineoutput_matrix');

%Produce scattered interpolants for vehicle data

%%
%OLDMODEL
% [MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
% M1_Grid = reshape(communicator(:,12),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% pres_Grid = reshape(communicator(:,13),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% temp_Grid = reshape(communicator(:,14),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% scattered.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
% scattered.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
% scattered.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');
%%
%NEWMODEL
shockdata = dlmread('ShockMat');
[MList,AOAList] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
scattered.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
scattered.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
scattered.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');
%%


scattered.data = dlmread('RESTM12DATA.txt');  
data = scattered.data;
% IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,6));
IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
M_englist = unique(sort(data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = floor(M_englist(1)):0.1:ceil(M_englist(end)); % enlarge spread, this is not necessary if you have a lot of engine data
% M_eng_interp = unique(sort(data(:,1)));

T_englist = unique(sort(data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = floor(T_englist(1)):10:ceil(T_englist(end)); 
% T_eng_interp = unique(sort(data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);
grid.Isp_eng = IspScattered(grid.Mgrid_eng,grid.T_eng);
scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','linear');
% scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'linear','linear');
% scattered.equivalence = scatteredInterpolant(data(:,1),data(:,2),data(:,7));
scattered.equivalence = scatteredInterpolant(data(:,1),data(:,2),data(:,4));


% Call LiftForceInterp
% This produces scattered interpolants which can calculate the vehicle
% settings required for trim at any flight conditions
[scattered.AoA,scattered.flapdeflection,scattered.drag,scattered.flap_pm,liftarray] = LiftForceInterp(communicator,communicator_trim,const,Atmosphere, scattered, SPARTAN_SCALE);
scattered.flap_def = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
scattered.flap_D = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
scattered.pm = scatteredInterpolant(communicator(:,1),communicator(:,2),communicator(:,11));
%==========================================================================
%   Import Payload Data 
%==========================================================================


% Import third stage data as arrays. the third stage data should be two
% files, named thirdstage.dat and thirdstagefine.dat, arranged in four
% columns: Altitude (m) , Trajectory angle (rad) , velocity (m/s) , payload-to-orbit (kg)

% Two payload arrays are used: course and fine, purely to minimise computational time 

% The PS routine must be able to search over a relatively large solution
% space for all primal variables end states, so there must be a
% payload-to-orbit solution at every possible end state.

ThirdStageData = dlmread('thirdstage.dat'); %Import Third Stage Data Raw 
ThirdStageData = sortrows(ThirdStageData);

% test of using griddedinterpolant
PayloadData = permute(reshape(ThirdStageData(:,6),[length(unique(ThirdStageData(:,5))),length(unique(ThirdStageData(:,4))),length(unique(ThirdStageData(:,3))),length(unique(ThirdStageData(:,2))),length(unique(ThirdStageData(:,1)))]),[5 4 3 2 1]);
[phiGrid,zetaGrid,VGrid,thetaGrid,vGrid] = ndgrid(unique(ThirdStageData(:,1)),unique(ThirdStageData(:,2)),unique(ThirdStageData(:,3)),unique(ThirdStageData(:,4)),unique(ThirdStageData(:,5)));
global PayloadGrid
PayloadGrid = griddedInterpolant(phiGrid,zetaGrid,VGrid,thetaGrid,vGrid,PayloadData,'spline','linear');

%=============================================== 

V0 = 20000.;
Vf = 50000.; %

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

VL = V0;
VU = 1.0*Vf; 

vL = 1500;
vU = 3100; % This limit must not cause the drag force to exceed the potential thrust of the vehicle by a large amount, otherwise DIDO will not solve

if const == 1  || const == 12 || const == 13 || const == 14
thetaL = -0.1;
else
thetaL = -0.1; %  NEED TO WATCH THAT THIS IS NOT OVERCONSTRAINING (ie. scramjet needs to be able to fly the optimal trajectory within these limits)
end

if const == 1  || const == 12 || const == 13 || const == 14
% thetaU = 0.1; % 
thetaU = 0.05; % 
else
thetaU = 0.1;  
end

mfuelL = 0;
mfuelU = 994*SPARTAN_SCALE; % 

% Define bounds of primals and controls ---------------------------------

%these must include every logical solution. For fuel mass I have found that
%bounds of exactly upper and lower fuel values over-constrain. I have
%modified the bounds accordingly


global scale

scale.V = 1;
scale.v = 1;
scale.theta = 1;
scale.thetadot = 1;
scale.m = 1;

if const == 1  || const == 12 || const == 14 || const == 13
bounds.lower.states = [VL/scale.V ; vL/scale.v; 0.1*thetaL/scale.theta; mfuelL/scale.m; -0.001/scale.thetadot];
bounds.upper.states = [VU/scale.V ; vU/scale.v; thetaU/scale.theta; (mfuelU+1)/scale.m; 0.002/scale.thetadot];

end

if const == 3 || const == 31
bounds.lower.states = [VL/scale.V ; vL/scale.v; thetaL/scale.theta; (mfuelL-3000)/scale.m; -0.001/scale.thetadot];
bounds.upper.states = [VU/scale.V ; vU/scale.v; thetaU/scale.theta; mfuelU/scale.m; 0.002/scale.thetadot];
end

% control bounds
if const == 1 || const == 12 || const == 13 || const == 14
omegadotL = -0.0001;
omegadotU = 0.0001;
% omegadotL = -0.00005;
% omegadotU = 0.00005;
else
% omegadotL = -0.001;
% omegadotU = 0.001;
% omegadotL = -0.0001;
% omegadotU = 0.0001;
omegadotL = -0.0001;
omegadotU = 0.0001;
end
bounds.lower.controls = [omegadotL/scale.thetadot];
bounds.upper.controls = [omegadotU/scale.thetadot]; 

%------------------
% bound the horizon
%------------------
% time bounds, this is unscaled

t0	    = 0;
tfMax 	= 450;   %  max tf; DO NOT set to Inf even for time-free problems % remember to set higher than Vmax bounds min time

% bounds.lower.time 	= [t0; t0];	
bounds.lower.time 	= [t0; 100];	
bounds.upper.time	= [t0; tfMax];

%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------

% v0 = 1764; 
v0 = 1524; 
% v0 = 1550; 
vf = 2839.51;

% See events file for definition of events function
if const == 1 || const == 12 || const == 14 || const == 13
    bounds.lower.events = [v0/scale.v; mfuelU/scale.m; mfuelL/scale.m]; % 
end


if const == 3 || const == 31
bounds.lower.events = [v0/scale.v; mfuelU/scale.m; mfuelL/scale.m];
end

bounds.upper.events = bounds.lower.events;      % equality event function bounds

%============================================
% Define the problem using DIDO expresssions:
%============================================
TwoStage2d.cost 		= 'SecondStageCost';
TwoStage2d.dynamics	    = 'SecondStageDynamics';
TwoStage2d.events		= 'SecondStageEvents';	

TwoStage2d.bounds       = bounds;


% Node Definition ====================================================
% the number of nodes is extremely important to the problem solution.
% usually between 50-150 works. The best node no. must be found using trial and error approach, but usually
% working down from 100 works well. 

if const == 3 || const == 31
% algorithm.nodes		= [60]; 
algorithm.nodes		= [80]; 
elseif const == 1
% algorithm.nodes		= [82];
algorithm.nodes		= [60]; 
elseif const == 12 
% algorithm.nodes		= [80];
algorithm.nodes		= [80];
elseif const == 13
% algorithm.nodes		= [78];
algorithm.nodes		= [80];
% algorithm.nodes		= [110];
elseif const == 14
algorithm.nodes		= [80];
end

global nodes

nodes = algorithm.nodes;


%  Guess =================================================================
constq = dlmread('primalconstq.txt');

% tfGuess = tfMax; % this needs to be close to make sure solution stays within Out_Force bounds

if const == 1

% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,34000 ];


guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,33000 ];


% guess.states(1,:) = constq(1,:);
elseif const == 12
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2)-100 ,34900]; %55kPa limited
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2)-100 ,35600];
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2)-100 ,32000];
elseif const == 13
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2)+100 ,34500];%45kPa limited
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2)-99 ,35800];%45kPa limited
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2)-100 ,34500];%45kPa limited
elseif const == 14
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,33500]; %High Drag

else
% guess.states(1,:) = [0 ,Vf]/scale.V; % for constant 50kPa
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/2860^2)+100];

guess.states(1,:) =[interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,32000 ]/scale.V; %50kpa limited
end

% guess.states(2,:) = [v0, vf]/scale.v; %v for normal use
guess.states(2,:) = [v0, 2700]/scale.v; %v for normal use


% guess.states(2,:) = constq(2,:);
if const ==3 || const == 31 
% guess.states(3,:) = [0,0]/scale.theta;
guess.states(3,:) = [0,0]/scale.theta;
% elseif const == 13
% guess.states(3,:) = [deg2rad(1.0),thetaU]/scale.theta;
% guess.states(3,:) = [deg2rad(0.5),0.04]/scale.theta;
else
% guess.states(3,:) = [0,0.05]/scale.theta;  
guess.states(3,:) = [deg2rad(1.3),0.05]/scale.theta;  

% guess.states(3,:) = constq(3,:);
end 

guess.states(4,:) = [mfuelU, 0]/scale.m;
% guess.states(4,:) = constq(4,:);

guess.states(5,:) = [0,0];
% guess.states(5,:) = zeros(1,length(constq(1,:)));

guess.controls(1,:)    = [0,0]; 
% guess.controls(1,:)    = zeros(1,length(constq(1,:)));

% guess.time        = [t0 ,tfMax];
% guess.time        = constq(7,:);
guess.time        = [t0 ,230];
% Tell DIDO the guess
%========================
algorithm.guess = guess;
% %========================
% algorithm.mode = 'accurate';
%=====================================================================================


%Start Plot
figure(10)
plot(linspace(guess.time(1),guess.time(2),algorithm.nodes),linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes),'k')
iterative_V(end+1,:) = linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes);
iterative_t(end+1,:) = linspace(guess.time(1),guess.time(2),algorithm.nodes);

filename = 'testnew51.gif';
frame = getframe(10);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',0);




% Call dido
% =====================================================================
tStart= cputime;    % start CPU clock 
[cost, primal, dual] = dido(TwoStage2d, algorithm);


% runTime = (cputime-tStart)
% runTime/60

EndTime = datestr(now,30)


V = primal.states(1,:)*scale.V;

v = primal.states(2,:)*scale.v;

t = primal.nodes;

theta = primal.states(3,:)*scale.theta;

mfuel = primal.states(4,:)*scale.m;

thetadot = primal.states(5,:)*scale.thetadot;

omegadot = primal.controls(1,:)*scale.theta;

% ===================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Stage
% Optimise third stage trajectory from end point for accuracy 
global zeta
global phi

cd('../ThirdStage')
ThirdStagePayloadMass = ThirdStageOptm(V(end),theta(end),v(end), phi(end),zeta(end))
cd('../SecondStage')


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
global Alpha
global a

Thrust = Thrust./cos(deg2rad(Alpha));

dt = t(2:end)-t(1:end-1); % Time change between each node pt
FuelUsed = zeros(1,nodes-1);
FuelUsed(1) = dt(1)*Fueldt(1);
for i = 2:nodes-1
    FuelUsed(i) = dt(i).*Fueldt(i) + FuelUsed(i-1);
end


% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(t(i+1) - t(i))*cos(theta(i)) + H(i);
end

Separation_LD = lift(end)/Fd(end)

figure(1)

subplot(5,5,[1,10])
hold on
plot(H, V)
plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')

dim = [.7 .52 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1000 - mfuel(end)) ' kg']},'FitBoxToText','on');  


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
plot(t, rad2deg(theta))
plot(t(algorithm.nodes(1)), rad2deg(theta(algorithm.nodes(1))), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory Angle (Deg)')



subplot(5,5,15)
plot(t, Fd)
title('Drag Force')

subplot(5,5,16)
hold on
plot(t, mfuel + 8755.1 - 994)
title('Vehicle Mass (kg)')



subplot(5,5,17)
plot(t, Thrust)
title('Thrust (N)')

Isp = Thrust./Fueldt./9.81;
IspNet = (Thrust-Fd)./Fueldt./9.81;

subplot(5,5,18)
plot(t, Isp)
title('Isp')

subplot(5,5,19)
plot(t, IspNet)
title('Net Isp')

subplot(5,5,20)
plot(t, flapdeflection)
title('Flap Deflection (deg)')

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

subplot(5,5,24)
hold on
plot(t, rad2deg(thetadot))
title('Trajectory Angle Change Rate (Deg/s)')

subplot(5,5,25)
hold on
plot(t, rad2deg(omegadot))
title('Omegadot Control (Deg/s2)')


dim = [.8 .0 .2 .2];
annotation('textbox',dim,'string',{['Third Stage Thrust: ', num2str(50), ' kN'],['Third Stage Starting Mass: ' num2str(2850) ' kg'],['Third Stage Isp: ' num2str(350) ' s']},'FitBoxToText','on');  

figure(2)
subplot(2,6,[1,6])
hold on
plot(H/1000, V/1000,'Color','k')

title('Trajectory')
xlabel('Horizontal Position (km)')
ylabel('Vertical Position (km)')

for i = 1:floor(t(end)/30)
    [j,k] = min(abs(t-30*i));
    str = strcat(num2str(round(t(k))), 's');
    text(H(k)/1000,V(k)/1000,str,'VerticalAlignment','top', 'FontSize', 10);
    
    plot(H(k)/1000, V(k)/1000, '+', 'MarkerSize', 10, 'MarkerEdgeColor','k')
end

plot(H(end)/1000, V(end)/1000, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','k')

text(H(end)/1000,V(end)/1000,'Third Stage Transition Point','VerticalAlignment','top', 'FontSize', 10);

dim = [.65 .45 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(994 - mfuel(end)) ' kg']},'FitBoxToText','on');  

thirdstageexample_H = [0+H(end) (H(end)-H(end - 1))+H(end) 20*(H(end)-H(end - 1))+H(end) 40*(H(end)-H(end - 1))+H(end) 60*(H(end)-H(end - 1))+H(end) 80*(H(end)-H(end - 1))+H(end)]/1000; %makes a small sample portion of an arbitrary third stage trajectory for example
thirdstageexample_V = [0+V(end) (V(end)-V(end - 1))+V(end) 20*((V(end)-V(end -1)))+V(end) 40*((V(end)-V(end -1)))+V(end) 60*((V(end)-V(end -1)))+V(end) 80*((V(end)-V(end -1)))+V(end)]/1000;
plot(thirdstageexample_H, thirdstageexample_V, 'LineStyle', '--','Color','k');

hold on
subplot(2,6,[7,9])
xlabel('time (s)')

hold on
ax1 = gca; % current axes
xlim([min(t) max(t)]);

line(t, rad2deg(theta),'Parent',ax1,'Color','k', 'LineStyle','-')

line(t, M,'Parent',ax1,'Color','k', 'LineStyle','--')

line(t, v./(10^3),'Parent',ax1,'Color','k', 'LineStyle','-.')

line(t, q./(10^4),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% line(t, heating_rate./(10^5),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% line(t, Q./(10^7),'Parent',ax1,'Color','k', 'LineStyle','-', 'lineWidth', 2.0)

% legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)',  'Q (Mj x 10)')
h = legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
rect1 = [0.12, 0.35, .25, .25];
set(h, 'Position', rect1)


subplot(2,6,[10,12])
xlabel('time (s)')
ax2 = gca;
xlim([min(t) max(t)]);
line(t, Alpha,'Parent',ax2,'Color','k', 'LineStyle','-')

line(t, flapdeflection,'Parent',ax2,'Color','k', 'LineStyle','--')


line(t, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')


line(t, IspNet./(10^2),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

figure(3)

subplot(2,5,[1,5]);

line(t, dual.dynamics(1,:),'Color','k', 'LineStyle','-');
line(t, dual.dynamics(2,:),'Color','k', 'LineStyle','--');
line(t, dual.dynamics(3,:),'Color','k', 'LineStyle','-.');
line(t, dual.dynamics(4,:),'Color','k', 'LineStyle',':');
line(t, dual.dynamics(5,:),'Color','k', 'LineStyle','-','LineWidth',2);
title('costates')
xlabel('time');
ylabel('Costates');
% axis([0,t(end),-1,1])
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4','\lambda_5');

subplot(2,5,[6,10])
Hamiltonian = dual.Hamiltonian(1,:);
plot(t,Hamiltonian,'Color','k');
axis([0,t(end),-1,1])
title('Hamiltonian')


% save results
dlmwrite('primal.txt', [primal.states;primal.controls;primal.nodes;q;IspNet;Alpha]);
dlmwrite('payload.txt', ThirdStagePayloadMass);
dlmwrite('dual.txt', [dual.dynamics;dual.Hamiltonian]);

copyfile('primal.txt',sprintf('../ArchivedResults/primal_%s.txt',Timestamp))
copyfile('dual.txt',sprintf('../ArchivedResults/dual_%s.txt',Timestamp))
copyfile('payload.txt',sprintf('../ArchivedResults/payload_%s.txt',Timestamp))

primal_old = primal;

ts = timeseries(Isp,t);
Mean_Isp = mean(ts)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% TESTING AND VALIDATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% If these are valid then solution is a KKT point

%COMPLEMENTARY CONDITIONS
% These should be zero if the state or control is within set bounds
% <=0 if at min bound, >=0 if at max bound

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
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This simply tests that the system dynamics hold, as the
% Pseudospectral method may not converge to a realistic
% solution



theta_F = cumtrapz(t,thetadot)+ theta(1);

thetadot_F = cumtrapz(t,omegadot) + thetadot(1);


v_F = cumtrapz(t,a);
v_F = v_F + v(1);

V_F = cumtrapz(t,v_F.*sin(theta_F));
V_F = V_F + V(1);

mfuel_F = cumtrapz(t,-Fueldt);
mfuel_F = mfuel_F + mfuel(1);

figure(6)

subplot(5,1,1)
plot(t,theta_F,t,theta);
title('Forward Simulation Comparison');
% note this is just a trapezoidal rule check, may not be exactly accurate
subplot(5,1,2)
plot(t,v_F,t,v);
subplot(5,1,3)
plot(t,V_F,t,V);
subplot(5,1,4)
plot(t,mfuel_F,t,mfuel);
subplot(5,1,4)
plot(t,thetadot_F,t,thetadot);

% Compute difference with CADAC for constant dynamic pressure path
if const == 3
    CADAC_DATA = dlmread('TRAJ.ASC');
    CADAC_Alpha = interp1(CADAC_DATA(:,1),CADAC_DATA(:,4),linspace(0,CADAC_DATA(end,1),nodes));
    CADAC_V = interp1(CADAC_DATA(:,1),CADAC_DATA(:,11),linspace(0,CADAC_DATA(end,1),nodes));
    MeanError_V = sum((CADAC_V - V)./V)/nodes
    MeanError_Alpha = sum((CADAC_Alpha - Alpha)./Alpha)/nodes
end



if PayloadGrid(phi(end),zeta(end),V(end)+10,theta(end),v(end)) - PayloadGrid(phi(end),zeta(end),V(end),theta(end),v(end)) < 0
    disp('Check Third Stage Payload Matrix, May Have Found False Maxima')
end


forward0 = [V(1),phi(1),theta(1),v(1),zeta(1),8755.1];


[f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);

figure(12)
hold on
plot(t,f_y(:,1));
plot(t,V);

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





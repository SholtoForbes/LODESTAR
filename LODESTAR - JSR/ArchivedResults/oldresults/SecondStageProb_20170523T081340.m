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

%%
% Copy the current setting to archive
% This saves the entire problem file every time the program is run. 

Timestamp = datestr(now,30)
copyfile('SecondStageProb.m',sprintf('../ArchivedResults/SecondStageProb_%s.m',Timestamp))
copyfile('SecondStageCost.m',sprintf('../ArchivedResults/SecondStageCost_%s.m',Timestamp))

%%
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
const = 12

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



%% Shock Data
% Import conical shock data and create interpolation splines 
shockdata = dlmread('ShockMat');
[MList,AOAList] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
scattered.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
scattered.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
scattered.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');
%% Engine Data
% Import engine data
scattered.data = dlmread('RESTM12DATA.txt');  
data = scattered.data;

newdata = [];
j=1;

% this sets the equivalence ratio interpolation region! VERY IMPORTANT
% The interpolators have trouble with equivalence ratio because its equal
% to 1 over a certain Mach no. (causes error in interpolator.

% This makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)
for i = 1: length(data(:,1))
    if data(i,1) < 5.
        newdata(j,:) = data(i,:);
        j=j+1;
    end
end

 p=polyfitn([data(:,1),data(:,2)],data(:,3),3);

M_englist = unique(sort(data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(data(:,1)));

T_englist = unique(sort(data(:,2))); % create unique list of angle of attack numbers from engine data 
T_eng_interp = unique(sort(data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

for i = 1:30
    for j= 1:30
grid.Isp_eng(i,j) = polyvaln(p,[grid.Mgrid_eng(i,j) grid.T_eng(i,j)]);
    end
end

scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');

scattered.equivalence = scatteredInterpolant(newdata(:,1),newdata(:,2),newdata(:,4), 'linear');
grid.eq_eng = scattered.equivalence(grid.Mgrid_eng,grid.T_eng);
scattered.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');




% Call LiftForceInterp
% This produces scattered interpolants which can calculate the vehicle
% settings required for trim at any flight conditions
% [scattered.AoA,scattered.flapdeflection,scattered.drag,scattered.flap_pm,liftarray] = LiftForceInterp(communicator,communicator_trim,const,Atmosphere, scattered, SPARTAN_SCALE);
liftarray = dlmread('liftarray');
[vList,altList,liftList] = ndgrid(unique(liftarray(:,1)),unique(liftarray(:,2)),unique(liftarray(:,3)));

AoA_Grid = permute(reshape(liftarray(:,4),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
scattered.AoA = griddedInterpolant(vList,altList,liftList,AoA_Grid,'spline','linear');

flapdeflection_Grid = permute(reshape(liftarray(:,5),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
scattered.flapdeflection = griddedInterpolant(vList,altList,liftList,flapdeflection_Grid,'spline','linear');

Drag_Grid = permute(reshape(liftarray(:,6),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
scattered.drag = griddedInterpolant(vList,altList,liftList,Drag_Grid,'spline','linear');

Flap_pitchingmoment_Grid = permute(reshape(liftarray(:,7),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
scattered.flap_pm = griddedInterpolant(vList,altList,liftList,Flap_pitchingmoment_Grid,'spline','linear');

scattered.flap_def = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
scattered.flap_D = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
scattered.pm = scatteredInterpolant(communicator(:,1),communicator(:,2),communicator(:,11));
clear liftarray

% scattered.drag_forward = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,4),liftarray(:,6));
% scattered.lift_forward = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,4),liftarray(:,3));
%==========================================================================
%   Import Payload Data 
%==========================================================================


% Import third stage data as arrays. the third stage data should be in thirdstage.dat
% columns: Altitude (m) , Trajectory angle (rad) , velocity (m/s) , payload-to-orbit (kg)

% The PS routine must be able to search over a relatively large solution
% space for all primal variables end states, so there must be a
% payload-to-orbit solution at every possible end state.

ThirdStageData = dlmread('thirdstage.dat'); %Import Third Stage Data Raw 
ThirdStageData = sortrows(ThirdStageData);



% This uses phi and zeta
% PayloadData = permute(reshape(ThirdStageData(:,6),[length(unique(ThirdStageData(:,5))),length(unique(ThirdStageData(:,4))),length(unique(ThirdStageData(:,3))),length(unique(ThirdStageData(:,2))),length(unique(ThirdStageData(:,1)))]),[5 4 3 2 1]);
% [phiGrid,zetaGrid,VGrid,thetaGrid,vGrid] = ndgrid(unique(ThirdStageData(:,1)),unique(ThirdStageData(:,2)),unique(ThirdStageData(:,3)),unique(ThirdStageData(:,4)),unique(ThirdStageData(:,5)));
% global PayloadGrid
% PayloadGrid = griddedInterpolant(phiGrid,zetaGrid,VGrid,thetaGrid,vGrid,PayloadData,'spline','linear');

%%Permutes payload matrix into interpolatable form, no phi, zeta

% PayloadData = permute(reshape(ThirdStageData(:,6),[length(unique(ThirdStageData(:,5))),length(unique(ThirdStageData(:,4))),length(unique(ThirdStageData(:,3)))]),[3 2 1]);
% [VGrid,thetaGrid,vGrid] = ndgrid(unique(ThirdStageData(:,3)),unique(ThirdStageData(:,4)),unique(ThirdStageData(:,5)));
% global PayloadGrid
% PayloadGrid = griddedInterpolant(VGrid,thetaGrid,vGrid,PayloadData,'spline','linear');


%THIS INTERPOLATEs FOR ANY MISSING THIRD STAGE POINTS, BE CAREFUL WITH THIS
[VGrid,thetaGrid,vGrid] = ndgrid(33000:1000:38000,[0 0.0125 0.025 0.0375 0.05],[2875:25:2950]);

PayloadDataInterp = scatteredInterpolant(ThirdStageData(:,3),ThirdStageData(:,4),ThirdStageData(:,5),ThirdStageData(:,6));
PayloadData = PayloadDataInterp(VGrid,thetaGrid,vGrid);
global PayloadGrid
PayloadGrid = griddedInterpolant(VGrid,thetaGrid,vGrid,PayloadData,'spline','linear');


% First Stage Array
FirstStageData = dlmread('FirstStageDat.txt');
scattered.FirstStagev = scatteredInterpolant(FirstStageData(:,2),FirstStageData(:,3),FirstStageData(:,4));

scattered.FirstStageData = FirstStageData; 

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

VL = 20000;
VU = 50000; 

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
% thetaU = 0.5; 
% thetaU = deg2rad(3); 
else
thetaU = 0.1;  
end

mfuelL = 0;
% mfuelU = 994*SPARTAN_SCALE; % 
 mfuelU = 1562*SPARTAN_SCALE; % 
% mfuelU = 1200*SPARTAN_SCALE; % 
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
bounds.lower.states = [VL/scale.V ; vL/scale.v; 0.1*thetaL/scale.theta; mfuelL/scale.m; -0.001/scale.thetadot; 1];
bounds.upper.states = [VU/scale.V ; vU/scale.v; thetaU/scale.theta; (mfuelU+1)/scale.m; 0.002/scale.thetadot; 2];

end

if const == 3 || const == 31
bounds.lower.states = [VL/scale.V ; vL/scale.v; thetaL/scale.theta; (mfuelL-3000)/scale.m; -0.001/scale.thetadot; 1];
bounds.upper.states = [VU/scale.V ; vU/scale.v; thetaU/scale.theta; mfuelU/scale.m; 0.002/scale.thetadot; 2];
end

% control bounds
if const == 1 || const == 12 || const == 13 || const == 14
% omegadotL = -0.0001;
% omegadotU = 0.0001;

omegadotL = -0.00005;
omegadotU = 0.00005;
% 
% omegadotL = -0.00005;
% omegadotU = 0.0001;
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
% tfMax 	= 450;   %  max tf; DO NOT set to Inf even for time-free problems % remember to set higher than Vmax bounds min time
 tfMax 	= 800;
% bounds.lower.time 	= [t0; t0];	
bounds.lower.time 	= [t0; 100];	
bounds.upper.time	= [t0; tfMax];

%-------------------------------------------
% Set up the bounds on the endpoint function
%-------------------------------------------

% v0 = 1764; 
% if const == 1 || const == 14 || const == 3
% v0 = 1524; 
% elseif const == 12
%     v0 = 1524;
% elseif const == 13
%     v0 = 1521;
% end


v0 = 1520; 
% v0 = 1490; % gives Mach 5.0018 for all cases
% v0 = 1600
vf = 2839.51;


% if const == 1 || const == 14
%     V0 = 24460;
%     
% elseif const == 12
%     V0 = 23860;
% 
% elseif const == 13
%     V0 = 25090;
% end


%% Define Events
% This defines set values along the trajectory.
% These correspond to the values in SecondStageEvents.m

% if const == 1 || const == 14 || const == 3
% bounds.lower.events = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2);v0/scale.v; mfuelU/scale.m; mfuelL/scale.m; 1.69];
% elseif const == 12
%     bounds.lower.events = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2);v0/scale.v; mfuelU/scale.m; mfuelL/scale.m; 1.69];
% elseif const == 13
%     bounds.lower.events = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2);v0/scale.v; mfuelU/scale.m; mfuelL/scale.m; 1.69];
% end


bounds.lower.events = [v0/scale.v; mfuelU/scale.m; mfuelL/scale.m; 1.69];


bounds.upper.events = bounds.lower.events;      % equality event function bounds


%% Define Path Constraints
% This limits the dynamic pressure.
if const == 1 || const == 14
%     bounds.lower.path = [0 ;0;0];
% bounds.upper.path = [50000 ;0;0];

    bounds.lower.path = 0;
bounds.upper.path = 50000;
elseif const == 12
%     bounds.lower.path = [0;0];
% bounds.upper.path = [55000;0];
    bounds.lower.path = 0;
bounds.upper.path = 55000;
elseif const == 13
%     bounds.lower.path = [0;0];
% bounds.upper.path = [45000;0];
    bounds.lower.path = 0;
bounds.upper.path = 45000;
elseif const ==3
%         bounds.lower.path = 0;
% bounds.upper.path = 0;
%         bounds.lower.path = [];
% bounds.upper.path = [];
end

%% 
%============================================
% Define the problem using DIDO expresssions:
%============================================
% Call the files whih DIDO uses
TwoStage2d.cost 		= 'SecondStageCost';
TwoStage2d.dynamics	    = 'SecondStageDynamics';
TwoStage2d.events		= 'SecondStageEvents';	
if const  == 1 || const  == 12 || const  == 13 || const  == 14 
TwoStage2d.path		= 'SecondStagePath';
end
TwoStage2d.bounds       = bounds;


%% Node Definition ====================================================
% The number of nodes is important to the problem solution. 
% While any number should theoretically work, in practice the choice of
% node number can have a large effect on results.
% Usually between 70-90 works. The best node no. must be found using trial and error approach, but usually
% working down from 90 works well. 

if const == 3 || const == 31
% algorithm.nodes		= [60]; 
algorithm.nodes		= [90]; 
elseif const == 1
algorithm.nodes		= [92]; 

% algorithm.nodes		= [89]; %works
elseif const == 12 
algorithm.nodes		= [89];
% algorithm.nodes		= [90]; % nearly works
elseif const == 13
% algorithm.nodes		= [78];
algorithm.nodes		= [90];
% algorithm.nodes		= [110];
elseif const == 14
algorithm.nodes		= [90];
end

global nodes

nodes = algorithm.nodes;


%%  Guess =================================================================
constq = dlmread('primalconstq.txt');

% tfGuess = tfMax; % this needs to be close to make sure solution stays within Out_Force bounds

% Altitude Guess. This is the most important guess, and should be close to
% the expected end solution. It is good for this end guess to be lower than
% expected.
if const == 1
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)+100 ,34500 ];
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,35000 ]; % this works ok
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-1000 ,35500 ];% works well
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-1000 ,36100 ];
elseif const == 12
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2)+100 ,35000];
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000/v0^2)-100 ,36000];
elseif const == 13
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2)+100 ,34500];%45kPa limited
 guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000/v0^2)-1000 ,36100];
elseif const == 14
% guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)+100 ,34000]; %High Drag
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100 ,35000]; 
else
% guess.states(1,:) =[interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2),33000 ]/scale.V; %50kpa limited
guess.states(1,:) =[interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2),32000 ]/scale.V; 
end

% Velocity Guess. This should be relatively close to the end solution,
% second most important guess. 
% guess.states(2,:) = [v0, 2950]/scale.v; 
if const == 1
%     guess.states(2,:) = [v0, 2920] % works
    guess.states(2,:) = [v0, 2940]
elseif const == 13
    guess.states(2,:) = [v0, 2900]
else
guess.states(2,:) = [v0, 2920]/scale.v; 
end
% Trajectoy angle guess. Sort of important, but easy to define. 
if const ==3 || const == 31 
guess.states(3,:) = [0,0]/scale.theta;
elseif const == 1 || const == 14
guess.states(3,:) = [0,0.05]/scale.theta;  
else
guess.states(3,:) = [0.05,0.05]
end 

% Mass guess. Simply the exact values at beginning and end (also constraints).
guess.states(4,:) = [mfuelU, 0]/scale.m;

% Trajectory angle change rate guess. Keep at 0.
guess.states(5,:) = [0,0];

% Heading angle guess. End defined. Start should be close to expected.
guess.states(6,:) = [1.682,1.699];

% Control guess. Keep at 0. Not important. 
guess.controls(1,:)    = [0,0]; 

% Time guess. This should be sort of close to expected (ie. not ridiculous)
if const == 1 || const == 14 || const == 3
guess.time        = [t0 ,350]; % works ok for 50kpa
else
guess.time        = [t0 ,400];
end
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

theta = primal.states(3,:)*scale.theta;

mfuel = primal.states(4,:)*scale.m;

thetadot = primal.states(5,:)*scale.thetadot;

zeta = primal.states(6,:);

omegadot = primal.controls(1,:)*scale.theta;

% ===================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Stage
% Optimise third stage trajectory from end point for accuracy 

global phi

cd('../ThirdStage')
[ThirdStagePayloadMass,ThirdStageControls,ThirdStageZeta,ThirdStagePhi,ThirdStageAlt,ThirdStagev,ThirdStaget,ThirdStageAlpha,ThirdStagem,ThirdStagegamma,ThirdStageq] = ThirdStageOptm(V(end),theta(end),v(end), phi(end),zeta(end));
ThirdStagePayloadMass
cd('../SecondStage')
% ThirdStagePayloadMass = 0;

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
global eq

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

figure(201)

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

figure(202)
subplot(2,6,[1,6])
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

text(H(end)/1000,V(end)/1000,'Third Stage Transition Point','VerticalAlignment','top', 'FontSize', 10);

dim = [.65 .45 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mfuel(1) - mfuel(end)) ' kg']},'FitBoxToText','on');  

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


% line(t, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
line(t, eq.*10,'Parent',ax2,'Color','k', 'LineStyle','-.')

line(t, IspNet./(10^2),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Equivalence Ratio x 10', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

figure(203)

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
dlmwrite('primal.txt', [primal.states;primal.controls;primal.nodes;q;IspNet;Alpha;M;eq;flapdeflection;phi]);
dlmwrite('payload.txt', ThirdStagePayloadMass);
dlmwrite('dual.txt', [dual.dynamics;dual.Hamiltonian]);
dlmwrite('ThirdStage.txt',[ThirdStageZeta;ThirdStagePhi;ThirdStageAlt;ThirdStagev;ThirdStaget;[ThirdStageAlpha 0];ThirdStagem;ThirdStagegamma;[ThirdStageq 0]]);

copyfile('primal.txt',sprintf('../ArchivedResults/primal_%s.txt',Timestamp))
copyfile('dual.txt',sprintf('../ArchivedResults/dual_%s.txt',Timestamp))
copyfile('payload.txt',sprintf('../ArchivedResults/payload_%s.txt',Timestamp))
copyfile('ThirdStage.txt',sprintf('../ArchivedResults/ThirdStage_%s.txt',Timestamp))
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

figure(205)

plot(t,dLHdu,t,mu_1,t,mu_2,t,mu_3,t,mu_4,t,mu_5,t,mu_u);
legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
title('Validation')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD INTEGRATION
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

figure(206)

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

t_diff = t - [0 t(1:end-1)];
if const == 3
    CADAC_DATA = dlmread('TRAJ.ASC');
    CADAC_Alpha = interp1(CADAC_DATA(:,2),CADAC_DATA(:,4),M(1:68)); % 1:68 gives mach numbers that align, may need to change this
    CADAC_V = interp1(CADAC_DATA(:,2),CADAC_DATA(:,11),M(1:68));
    MeanError_V = sum(abs(CADAC_V - V(1:68))./V(1:68).*t_diff(1:68))/t(end)
    MeanError_Alpha = sum(abs(CADAC_Alpha - Alpha(1:68))./Alpha(1:68).*t_diff(1:68))/t(end)
end

% if PayloadGrid(phi(end),zeta(end),V(end)+10,theta(end),v(end)) - PayloadGrid(phi(end),zeta(end),V(end),theta(end),v(end)) < 0
%     disp('Check Third Stage Payload Matrix, May Have Found False Maxima')
% end
if PayloadGrid(V(end)+10,theta(end),v(end)) - PayloadGrid(V(end),theta(end),v(end)) < 0
    disp('Check Third Stage Payload Matrix, Found Maxima')
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 

forward0 = [V(1),phi(1),theta(1),v(1),zeta(1),9.7725e+03];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered,AlphaInterp(t,lift,f_t),AlphaInterp(t,Fd,f_t),AlphaInterp(t,Thrust,f_t),AlphaInterp(t,flapdeflection,f_t)),t(1:end),forward0);

figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(t,V);


%% plot engine interpolation visualiser
T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), V); 
T1 = scattered.tempgridded(M,Alpha).*T0;
M1 = scattered.M1gridded(M, Alpha);

plotM = [min(M_englist):0.01:10];
plotT = [min(T_englist):1:600];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = scattered.eqGridded(gridM,gridT);
interpIsp = scattered.IspGridded(gridM,gridT);

figure(210)
hold on
contourf(gridM,gridT,interpeq);
scatter(data(:,1),data(:,2),30,data(:,4),'filled');
plot(M1,T1,'r');

error_Isp = scattered.IspGridded(data(:,1),data(:,2))-data(:,3);

figure(211)
hold on
contourf(gridM,gridT,interpIsp);
scatter(data(:,1),data(:,2),30,data(:,3),'filled')
plot(M1,T1,'r');

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);



% Run First Stage =========================================================
cd('../FirstStage')
[FirstStageStates] = FirstStageProblem(V(1),theta(1),phi(1),zeta(1),const);
cd('../SecondStage')
dlmwrite('FirstStage.txt', FirstStageStates);
copyfile('FirstStage.txt',sprintf('../ArchivedResults/FirstStage_%s.txt',Timestamp))
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






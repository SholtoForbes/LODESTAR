%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\..\thirdStage-GPOPS')
addpath('..\EngineData')
addpath('..\')
%% Atmosphere Data %%======================================================
Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data

auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data

auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data

auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 

auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 

%% Import Vehicle and trajectory Config Data %%============================
addpath('../')
run VehicleConfig.m
run TrajectoryConfig50kPa.m

auxdata.Stage3 = Stage3;
auxdata.Stage2 = Stage2;


%%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)

auxdata.A = 62.77; %m^2
%%
% Copy the current setting to archive
% This saves the entire problem file every time the program is run. 

% Timestamp = datestr(now,30)
% mkdir('../ArchivedResults', sprintf(Timestamp))
% copyfile('SecondStageProb.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',Timestamp))
% copyfile('SecondStageCost.m',sprintf('../ArchivedResults/%s/SecondStageCost.m',Timestamp))

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
% 32: Higher velocity


const = 1
auxdata.const = const;
%% Aerodynamic Data - Communicator %%======================================
% Take inputs of aerodynamic communicator matrices, these should be .txt files 
% This is used for forward simulation. 
% communicator = importdata('communicator.txt');
% communicator_trim = importdata('communicator_trim.txt');
% 
% auxdata.interp.flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
% auxdata.interp.flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
% auxdata.interp.flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,6));
% 
% [MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
% Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% 
% auxdata.interp.Cl_spline1 = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
% auxdata.interp.Cd_spline1 = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
% auxdata.interp.pitchingmoment_spline1 = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');

%% Aerodynamic Data 


addpath ..\CG14.5
% 

aero_EngineOff.fullFuel = importdata('SPARTANaero14.9122');
flapaero.fullFuel = importdata('SPARTANaeroFlaps14.9122');
aero_EngineOn.fullFuel = importdata('SPARTANaeroEngineOn14.9122');

aero_EngineOff.noFuel = importdata('SPARTANaero15.3515');
flapaero.noFuel = importdata('SPARTANaeroFlaps15.3515');
aero_EngineOn.noFuel = importdata('SPARTANaeroEngineOn15.3515');

[auxdata.interp.Cl_spline_EngineOff.fullFuel,auxdata.interp.Cd_spline_EngineOff.fullFuel,auxdata.interp.Cl_spline_EngineOn.fullFuel,auxdata.interp.Cd_spline_EngineOn.fullFuel,auxdata.interp.flap_spline_EngineOff.fullFuel,auxdata.interp.flap_spline_EngineOn.fullFuel] = AeroInt(aero_EngineOff.fullFuel,aero_EngineOn.fullFuel,flapaero.fullFuel);

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero_EngineOff.noFuel,aero_EngineOn.noFuel,flapaero.noFuel);


%% Conical Shock Data %%===================================================
% Import conical shock data and create interpolation splines 
shockdata = dlmread('ShockMat');
[MList,AOAList] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
auxdata.interp.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
auxdata.interp.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');


%% Engine Data %%==========================================================
% Import engine data
auxdata.interp.engine_data = dlmread('ENGINEDATA.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.interp.engine_data;

% Create uniform grid of Mach no. and temperature values. 
M_englist = unique(sort(engine_data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(engine_data(:,1)));

T_englist = unique(sort(engine_data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = unique(sort(engine_data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

% Set the equivalence ratio interpolation region %-------------------------
% VERY IMPORTANT

% The interpolators have trouble with equivalence ratio because its equal
% to 1 over a certain Mach no. (causes error in interpolator, as the
% interpolator will find values of equivalence ratio < 1 where they should
% not exist)

% This makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)

% the the maximum of this to around where equivalence ratio stops changing,
% and check the end results

eq_data = [];
j=1;
for i = 1: length(engine_data(:,1))
    if engine_data(i,1) < 5.
        eq_data(j,:) = engine_data(i,:);
        j=j+1;
    end
end

auxdata.interp.equivalence = scatteredInterpolant(eq_data(:,1),eq_data(:,2),eq_data(:,4), 'linear');
grid.eq_eng = auxdata.interp.equivalence(grid.Mgrid_eng,grid.T_eng);
auxdata.interp.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');

% Load the interpolated Isp data %-----------------------------------------

% gridIsp_eng is the spline interpolated data set created by
% engineint.m and engineinterpolator.exe

load gridIsp_eng
grid.Isp_eng = gridIsp_eng;

% gridIsp_eng may have sections at which the Isp is 0. The following finds
% these, and fills them in with linearly intepolated values.
Isp_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,3));

for i = 1:30 % must match engineint.m
    for j= 1:30
        % grid.Isp_eng(i,j) = polyvaln(p,[grid.Mgrid_eng(i,j) grid.T_eng(i,j)]);
        if any(grid.Isp_eng(i,j)) == false
            grid.Isp_eng(i,j) = Isp_interpolator(grid.Mgrid_eng(i,j), grid.T_eng(i,j));
        end
    end
end

auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');


%% Import Payload Data %%==================================================

% Import third stage data as arrays. the third stage data should be in thirdstage.dat
% columns: Altitude (m) , Trajectory angle (rad) , velocity (m/s) , payload-to-orbit (kg)

% The PS routine must be able to search over a relatively large solution
% space for all primal variables end states, so there must be a
% payload-to-orbit solution at every possible end state.

ThirdStageData = dlmread('thirdstage.dat'); %Import Third Stage Data Raw 
ThirdStageData = sortrows(ThirdStageData);

% Interpolate for Missing Third Stage Points %-----------------------------
% Be careful with this. 
[VGrid,gammaGrid,vGrid] = ndgrid(unique(ThirdStageData(:,3)),unique(ThirdStageData(:,4)),unique(ThirdStageData(:,5))); % must match the data in thirdstage.dat

PayloadDataInterp = scatteredInterpolant(ThirdStageData(:,3),ThirdStageData(:,4),ThirdStageData(:,5),ThirdStageData(:,6)); % interpolate for missing third stage points

PayloadData = PayloadDataInterp(VGrid,gammaGrid,vGrid);

auxdata.PayloadGrid = griddedInterpolant(VGrid,gammaGrid,vGrid,PayloadData,'spline','linear');

%% Import Bounds %%========================================================
lonMin = -pi;         lonMax = -lonMin;
latMin = -70*pi/180;  latMax = -latMin;
lat0 = -0.264;
lon0 = deg2rad(145);
aoaMin = deg2rad(0);  aoaMax = 9*pi/180;
% bankMin1 = -1*pi/180; bankMax1 =   50*pi/180;

% Primal Bounds
bounds.phase(1).state.lower = [Stage2.Bounds.Alt(1), lonMin, latMin, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin, Stage2.Bounds.mFuel(1)];
bounds.phase(1).state.upper = [Stage2.Bounds.Alt(2), lonMax, latMax, Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax, Stage2.Bounds.mFuel(2)];

% Initial States
bounds.phase(1).initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, Stage2.Initial.v, Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin, Stage2.Initial.mFuel] ;
bounds.phase(1).initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, Stage2.Initial.v, Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax, Stage2.Initial.mFuel];

% End States
bounds.phase(1).finalstate.lower = [Stage2.Bounds.Alt(1), lonMin, latMin, Stage2.Bounds.v(1), Stage2.End.gammaOpt(1), Stage2.End.Zeta, aoaMin, Stage2.End.mFuel];
bounds.phase(1).finalstate.upper = [Stage2.Bounds.Alt(2), lonMax, latMax, Stage2.Bounds.v(2), Stage2.End.gammaOpt(2), Stage2.End.Zeta, aoaMax, Stage2.Initial.mFuel];

% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-.1)];
bounds.phase(1).control.upper = [deg2rad(.1)];
% Time Bounds

bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = Stage2.Bounds.time(1);
bounds.phase(1).finaltime.upper = Stage2.Bounds.time(2);

%% Define Path Constraints
% This limits the dynamic pressure.
if const == 1 || const == 14 || const == 15
    bounds.phase(1).path.lower = [0];
    bounds.phase(1).path.upper = [50000];
elseif const == 12
    bounds.phase(1).path.lower = [0 ];
    bounds.phase(1).path.upper = [55000 ];
elseif const == 13
    bounds.phase(1).path.lower = [0 ];
    bounds.phase(1).path.upper = [45000];
elseif const ==3 || const == 32
        bounds.phase(1).path.lower = [0 ];
    bounds.phase(1).path.upper = [50010];
end


%%  Guess =================================================================

guess.phase(1).state(:,1)   = [22000;25000];
guess.phase(1).state(:,2)   = [0;0];
guess.phase(1).state(:,3)   = [-0.269;-0.13];
guess.phase(1).state(:,4)   = Stage2.Guess.v.';
guess.phase(1).state(:,5)   = Stage2.Guess.gamma.';
guess.phase(1).state(:,6)   = Stage2.Guess.zeta.';
guess.phase(1).state(:,7)   = [2*pi/180; 2*pi/180];
guess.phase(1).state(:,8) 	= [Stage2.Initial.mFuel, 200];

guess.phase(1).control      = [[0;0]];
guess.phase(1).time          = [0;650];

% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 



%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
% mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 5;
mesh.colpointsmin = 3;
mesh.colpointsmax = 50;
mesh.tolerance    = 1e-5;


%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @SecondStageContinuous;
setup.functions.endpoint             = @SecondStageEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.maxiterations = 1000;
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
% setup.scales.method                  = 'automatic-guessUpdate';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%


output = gpops2(setup);

%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================
% Assign the primal variables
alt = output.result.solution.phase(1).state(:,1).';
lon = output.result.solution.phase(1).state(:,2).';
lat = output.result.solution.phase(1).state(:,3).';
v = output.result.solution.phase(1).state(:,4).'; 
gamma = output.result.solution.phase(1).state(:,5).'; 
zeta = output.result.solution.phase(1).state(:,6).';
Alpha = output.result.solution.phase(1).state(:,7).';
% eta = output.result.solution.phase(1).state(:,8).';
mFuel = output.result.solution.phase(1).state(:,8).'; 

omegadot  = output.result.solution.phase(1).control.'; 


time = output.result.solution.phase(1).time.';

figure(201)
subplot(9,1,1)
hold on
plot(time,alt)
subplot(9,1,2)
hold on
plot(time,v)
subplot(9,1,3)
hold on
plot(time,lon)
subplot(9,1,4)
hold on
plot(time,lat)
subplot(9,1,5)
hold on
plot(time,v)
subplot(9,1,6)
hold on
plot(time,gamma)
subplot(9,1,7)




% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point

global phi

% cd('../ThirdStage')
% [ThirdStagePayloadMass,ThirdStageControls,ThirdStageZeta,ThirdStagePhi,ThirdStageAlt,ThirdStagev,ThirdStaget,ThirdStageAlpha,ThirdStagem,ThirdStagegamma,ThirdStageq] = ThirdStageOptm(V(end),gamma(end),v(end), phi(end),zeta(end), 1);
% ThirdStagePayloadMass
% cd('../SecondStage')
ThirdStagePayloadMass = 0;

[~,~,~,~,~,~, q, M, D, rho,L,Fueldt,T,Isp,flapdeflection] = VehicleModelAscent(gamma, alt, v,auxdata,zeta,lat,lon,Alpha,0,1,mFuel,mFuel(1),mFuel(end),1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt)

% eq = Engine.eq;
% Thrust = Engine.Thrust;
% Fueldt = Engine.Fueldt;

% Fd = Vehicle.Fd;
% Alpha = Vehicle.Alpha;
% lift = Vehicle.lift;
% flapdeflection = Vehicle.flapdeflection;
% 
% Thrust = Thrust./cos(deg2rad(Alpha)); % change thrust to account for total thrust, including portion that contributes to lift
% 
% dt = time(2:end)-time(1:end-1); % Time change between each node pt
% FuelUsed = zeros(1,nodes-1);
% FuelUsed(1) = dt(1)*Fueldt(1);
% for i = 2:nodes-1
%     FuelUsed(i) = dt(i).*Fueldt(i) + FuelUsed(i-1);
% end


% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
end

% Separation_LD = lift(end)/Fd(end)

figure(201)

subplot(5,5,[1,10])
hold on
plot(H, alt)
% plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')

dim = [.7 .52 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1000 - mFuel(end)) ' kg']},'FitBoxToText','on');  


subplot(5,5,11)
hold on
plot(time, v)

title('Velocity (m/s)')


subplot(5,5,12)
plot(time, M)
title('Mach no')

subplot(5,5,13)
plot(time, q)
title('Dynamic Pressure (pa)')

subplot(5,5,14)
hold on
plot(time, rad2deg(gamma))

title('Trajectory Angle (Deg)')



subplot(5,5,15)
plot(time, D)
title('Drag Force')

subplot(5,5,16)
hold on
plot(time, mFuel + 8755.1 - 994)
title('Vehicle Mass (kg)')



subplot(5,5,17)
plot(time, T)
title('Thrust (N)')

% Isp = Thrust./Fueldt./9.81;
IspNet = (T-D)./Fueldt./9.81;

subplot(5,5,18)
plot(time, Isp)
title('Isp')

subplot(5,5,19)
plot(time, IspNet)
title('Net Isp')

subplot(5,5,20)
plot(time, flapdeflection)
title('Flap Deflection (deg)')

subplot(5,5,21)
plot(time, rad2deg(Alpha))
title('Angle of Attack (deg)')

% subplot(5,5,22);
% plot(time, dual.dynamics);
% title('costates')
% xlabel('time');
% ylabel('costates');
% legend('\lambda_1', '\lambda_2', '\lambda_3');

% subplot(5,5,23)
% Hamiltonian = dual.Hamiltonian(1,:);
% plot(time,Hamiltonian);
% title('Hamiltonian')

% subplot(5,5,24)
% hold on
% plot(time, rad2deg(gammadot))
% title('Trajectory Angle Change Rate (Deg/s)')
% 
% subplot(5,5,25)
% hold on
% plot(time, rad2deg(omegadot))
% title('Omegadot Control (Deg/s2)')


dim = [.8 .0 .2 .2];
annotation('textbox',dim,'string',{['Third Stage Thrust: ', num2str(50), ' kN'],['Third Stage Starting Mass: ' num2str(2850) ' kg'],['Third Stage Isp: ' num2str(350) ' s']},'FitBoxToText','on');  

figure(202)
sp1 = subplot(2,6,[1,6]);
ax1 = gca; % current axes
hold on
plot(H/1000, alt/1000,'Color','k')

title('Trajectory')
xlabel('Earth Normal Distance Flown (km)')
ylabel('Vertical Position (km)')

for i = 1:floor(time(end)/30)
    [j,k] = min(abs(time-30*i));
    str = strcat(num2str(round(time(k))), 's');
    text(H(k)/1000,alt(k)/1000,str,'VerticalAlignment','top', 'FontSize', 10);
    
    plot(H(k)/1000, alt(k)/1000, '+', 'MarkerSize', 10, 'MarkerEdgeColor','k')
end

plot(H(end)/1000, alt(end)/1000, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','k')

text(H(end)/1000,alt(end)/1000,'Third Stage Transition Point','VerticalAlignment','top', 'FontSize', 10);

dim = [.65 .45 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mFuel(1) - mFuel(end)) ' kg']},'FitBoxToText','on');  

thirdstageexample_H = [0+H(end) (H(end)-H(end - 1))+H(end) 20*(H(end)-H(end - 1))+H(end) 40*(H(end)-H(end - 1))+H(end) 60*(H(end)-H(end - 1))+H(end) 80*(H(end)-H(end - 1))+H(end)]/1000; %makes a small sample portion of an arbitrary third stage trajectory for example
thirdstageexample_V = [0+alt(end) (alt(end)-alt(end - 1))+alt(end) 20*((alt(end)-alt(end -1)))+alt(end) 40*((alt(end)-alt(end -1)))+alt(end) 60*((alt(end)-alt(end -1)))+alt(end) 80*((alt(end)-alt(end -1)))+alt(end)]/1000;
plot(thirdstageexample_H, thirdstageexample_V, 'LineStyle', '--','Color','k');

hold on
sp2 = subplot(2,6,[7,9]);
xlabel('time (s)')

hold on
ax2 = gca; % current axes
xlim([min(time) max(time)]);

line(time, rad2deg(gamma),'Parent',ax2,'Color','k', 'LineStyle','-')

line(time, M,'Parent',ax2,'Color','k', 'LineStyle','--')

line(time, v./(10^3),'Parent',ax2,'Color','k', 'LineStyle','-.')

line(time, q./(10^4),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% line(time, heating_rate./(10^5),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% line(time, Q./(10^7),'Parent',ax1,'Color','k', 'LineStyle','-', 'lineWidth', 2.0)

% legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)',  'Q (Mj x 10)')
h = legend(ax2,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
rect1 = [0.12, 0.35, .25, .25];
set(h, 'Position', rect1)


sp3 = subplot(2,6,[10,12]);
xlabel('time (s)')
ax3 = gca;
xlim([min(time) max(time)]);
line(time, rad2deg(Alpha),'Parent',ax3,'Color','k', 'LineStyle','-')

line(time, flapdeflection,'Parent',ax3,'Color','k', 'LineStyle','--')


% line(time, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
% line(time, eq.*10,'Parent',ax3,'Color','k', 'LineStyle','-.')

% line(time, IspNet./(10^2),'Parent',ax3,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax3, 'AoA (degrees)','Flap Deflection (degrees)', 'Equivalence Ratio x 10', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

if auxdata.PayloadGrid(alt(end)+10,gamma(end),v(end)) - auxdata.PayloadGrid(alt(end),gamma(end),v(end)) < 0
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




forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),0,1,mFuel(1),mFuel(end)),time(1:end),forward0);


figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,gamma);

subplot(7,1,5)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time,v);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time,zeta);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time,mFuel);

%% plot engine interpolation visualiser
T0 = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2), alt); 
T1 = interp.tempgridded(M,Alpha).*T0;
M1 = interp.M1gridded(M, Alpha);

plotM = [min(M_englist):0.01:9.5];
plotT = [min(T_englist):1:550];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = interp.eqGridded(gridM,gridT);
interpIsp = interp.IspGridded(gridM,gridT);

figure(210)
hold on
contourf(gridM,gridT,interpeq);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'filled');
plot(M1,T1,'r');

error_Isp = interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);

figure(211)
hold on
contourf(gridM,gridT,interpIsp);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'filled')
plot(M1,T1,'r');

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);



% Run First Stage =========================================================
cd('../FirstStage')
[FirstStageStates] = FirstStageProblem(alt(1),gamma(1),phi(1),zeta(1),const);
cd('../SecondStage')
dlmwrite('FirstStage.txt', FirstStageStates);
copyfile('FirstStage.txt',sprintf('../ArchivedResults/%s/firststage_%s.txt',Timestamp,Timestamp))


%% Latitude Plot
figure(250)
plot(FirstStageStates(:,9))
plot(phi)
plot(ThirdStagePhi)
title('Latitude')

%% SAVE FIGS
saveas(figure(301),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ThirdStage.fig']);
saveas(figure(101),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'FirstStage.fig']);
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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc

%%
% =========================================================================
% SET RUN MODE
% =========================================================================
% Change mode to set the target of the simulation. Much of the problem
% definition changes with mode.

% mode = 1x: No end constraint, used for optimal trajectory calculation
% mode = 1: 50kPa limit, 12: 55 kPa limit, 13: 45 kPa limit, 14: 50kPa limit & 10% additional drag

mode = 1
% auxdata.mode = mode;

%% Misc Modifiers

auxdata.delta = deg2rad(0) % thrust vector angle test
auxdata.dragmod = 1. %drag increase test

%% Launch Point
lat0 = deg2rad(-12.164); % Equatorial Launch Australia Spaceport near Nhulunbuy
lon0 = deg2rad(136.755);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\..\thirdStage')
addpath('..\EngineData')
addpath('..\')

Timestamp = datestr(now,30)
mkdir('../ArchivedResults', strcat(Timestamp, 'mode', num2str(mode)))
copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',strcat(Timestamp,'mode',num2str(mode))))

%% Atmosphere Data %%======================================================
% Fetch atmospheric data and compute interpolation splines.

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

%% Third Stage Aerodynamic Data

auxdata.Aero3 = dlmread('aerocoeffs.txt');

Aero3 = auxdata.Aero3;
auxdata.interp.Drag_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,5));
% 
auxdata.interp.Lift_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,6));

auxdata.interp.CP_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,7));

auxdata.interp.CN_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,4));

auxdata.interp.Max_AoA_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,4),Aero3(:,2));


%% Conical Shock Data %%===================================================
% Import conical shock data and create interpolation splines 
shockdata = dlmread('ShockMat');
[MList_EngineOn,AOAList_EngineOn] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
auxdata.interp.M1gridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,M1_Grid,'spline','linear');
auxdata.interp.presgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,pres_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,temp_Grid,'spline','linear');

%% Equivalence Ratio %%==========================================================
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

%% Isp data %-----------------------------------------

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
            grid.Isp_eng(i,j) = Isp_interpolator(grid.Mgrid_eng(i,j), grid.T_eng(i,j)); % Linearly extrapolate for any data point which the Bivar spline could not solve
        end
    end
end

auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');



%% Aerodynamic Data

% Fetch aerodynamic data and compute interpolation splines.
% Calculate the flap deflection necessary for trim.
% Each set of aero corresponds to a different CG. 

addpath ..\ViscousAero

Viscousaero_EngineOn = importdata('VC3D_viscousCoefficients_ascent.dat');
Viscousaero_EngineOff = importdata('VC3D_viscousCoefficients_descent.dat');

% These aerodynamic datasets have been created in ClicCalcCGVar.m

T_L = -1.327; % Thrust location in z, average (m), measured from CREO

addpath ..\CG15.1255
% Full of fuel, with third stage

CG_z = (-0.1974*(4.9571e+03+1562) + 3300*0.547)/(4.9571e+03+1562+3300);

aero1.aero_EngineOff = importdata('SPARTANaero15.228');
aero1.flapaero = importdata('SPARTANaeroFlaps15.228');
aero1.aero_EngineOn = importdata('SPARTANaeroEngineOn15.228');
aero1.aero_Engine = importdata('SPARTANEngine15.228');
aero1.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero1.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.fullFuel,auxdata.interp.Cd_spline_EngineOff.fullFuel,auxdata.interp.Cl_spline_EngineOn.fullFuel,auxdata.interp.Cd_spline_EngineOn.fullFuel,auxdata.interp.flap_spline_EngineOff.fullFuel,auxdata.interp.flap_spline_EngineOn.fullFuel] = AeroInt(aero1,auxdata,T_L,CG_z);


% Cylindrical fuel tanks depleted, with third stage
CG_z = (-0.1974*(4.9571e+03+852) + 3300*0.547)/(4.9571e+03+852+3300);

aero1.aero_EngineOff = importdata('SPARTANaero15.1557');
aero1.flapaero = importdata('SPARTANaeroFlaps15.1557');
aero1.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1557');
% aero1.aero_Engine = importdata('SPARTANEngine15.1557');
aero1.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero1.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.cylTankEnd,auxdata.interp.Cd_spline_EngineOff.cylTankEnd,auxdata.interp.Cl_spline_EngineOn.cylTankEnd,auxdata.interp.Cd_spline_EngineOn.cylTankEnd,auxdata.interp.flap_spline_EngineOff.cylTankEnd,auxdata.interp.flap_spline_EngineOn.cylTankEnd] = AeroInt(aero1,auxdata,T_L,CG_z);



% End of acceleration, with third stage. 
% NOTE: it is assumed that the CG does not change
% due to fuel after the end of acceleration phase, so there will still be some fuel left when this is used. 

CG_z = (-0.2134*4.9571e+03+ 3300*0.547)/(4.9571e+03+3300);

aero2.aero_EngineOff = importdata('SPARTANaero15.727');
aero2.flapaero = importdata('SPARTANaeroFlaps15.727');
aero2.aero_EngineOn = importdata('SPARTANaeroEngineOn15.727');
aero2.aero_Engine = importdata('SPARTANEngine15.727');
aero2.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero2.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero2,auxdata,T_L,CG_z);

% Flyback, without third stage. Fuel variation not used for flyback.

CG_z = -0.2134; % calculated fom CREO

aero3.aero_EngineOff = importdata('SPARTANaero15.1255');
aero3.flapaero = importdata('SPARTANaeroFlaps15.1255');
aero3.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1255');
aero3.aero_Engine = importdata('SPARTANEngine15.1255');
aero3.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero3.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStage,auxdata.interp.Cd_spline_EngineOff.noThirdStage,auxdata.interp.Cl_spline_EngineOn.noThirdStage,auxdata.interp.Cd_spline_EngineOn.noThirdStage,auxdata.interp.flap_spline_EngineOff.noThirdStage,auxdata.interp.flap_spline_EngineOn.noThirdStage] = AeroInt(aero3,auxdata,T_L,CG_z);


%% Import Bounds %%========================================================

latMin2 = -0.5;  latMax2 = 0.5;

% lat0 = -0.264;
% lon0 = deg2rad(145);

aoaMin21 = 0;  aoaMax21 = 10*pi/180;
bankMin21 = -1*pi/180; 
% bankMax21 =   90*pi/180;
bankMax21 =   1*pi/180;
% Primal Bounds
bounds.phase.state.lower = [Stage2.Bounds.Alt(1), 2, -1, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Bounds.mFuel(1)];
bounds.phase.state.upper = [Stage2.Bounds.Alt(2), 3, 1, Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Bounds.mFuel(2)];

% Initial States
bounds.phase.initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, 1500, 0, deg2rad(90), aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
bounds.phase.initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, 1510, deg2rad(15), deg2rad(90), aoaMax21, bankMax21, Stage2.Initial.mFuel];
% bounds.phase.initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, 1500, Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
% bounds.phase.initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, 1500, deg2rad(15), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Initial.mFuel];

% End States
% End bounds are set slightly differently, to encourage an optimal solution
bounds.phase.finalstate.lower = [20000, 2, -1, 2800, 0, Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.End.mFuel];
bounds.phase.finalstate.upper = [45000, 3, 1, Stage2.Bounds.v(2), deg2rad(20), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.End.mFuel];

% bounds.phase.finalstate.lower = [34000, lonMin, latMin2, 2300, 0, Stage2.Bounds.zeta(1), aoaMin21, 0, Stage2.End.mFuel];
% bounds.phase.finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), 0, Stage2.Bounds.zeta(2), aoaMax21, 0, Stage2.Initial.mFuel];
%  disp('gamma set to 0')
 
% Control Bounds
bounds.phase.control.lower = [deg2rad(-.5), deg2rad(-1)];
bounds.phase.control.upper = [deg2rad(.5), deg2rad(1)];

% Time Bounds
bounds.phase.initialtime.lower = 0;
bounds.phase.initialtime.upper = Stage2.Bounds.time(2);
bounds.phase.finaltime.lower = Stage2.Bounds.time(1);
bounds.phase.finaltime.upper = Stage2.Bounds.time(2);


%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.
if mode == 1 || mode == 14 || mode == 15
    bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [70000];
elseif mode == 12
    bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [55000];
elseif mode == 13
    bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [45000];
elseif mode ==3 || mode == 32
        bounds.phase.path.lower = [50000];
    bounds.phase.path.upper = [51000];
end

bounds.phase.integral.lower = -10000000000;
bounds.phase.integral.upper = 10000000000;
%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
% guess.phase.state(:,1)   = [24000;35000];
guess.phase.state(:,2)   = [2.53;2.5368];
guess.phase.state(:,3)   = [-0.269;-0.10];
guess.phase.state(:,4)   = Stage2.Guess.v.';
guess.phase.state(:,5)   = Stage2.Guess.gamma.';
guess.phase.state(:,6)   = Stage2.Guess.zeta.';
guess.phase.state(:,7)   = [2*pi/180; 5*pi/180];
guess.phase.state(:,8)   = [deg2rad(10);deg2rad(10)];
guess.phase.state(:,9) 	= [Stage2.Initial.mFuel; 100];

guess.phase.control      = [[0;0],[0;0]];
guess.phase.time          = [0;650];

guess.phase.integral = 0

% Tie stages together
bounds.eventgroup.lower = [0];
bounds.eventgroup.upper = [1500]; 

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre'; % Default method does not perform h adaptions sometimes, and can be unstable. 
%  mesh.method       = 'hp-DarbyRao';
mesh.maxiterations = 3;
mesh.colpointsmin = 8;
mesh.colpointsmax = 50;
mesh.tolerance    = 1e-5;

% mesh.phase.fraction = 0.1*ones(1,10)
% 
% mesh.phase(3).fraction = 0.05*ones(1,20)
% 
% mesh.phase(4).fraction = 0.1*ones(1,10)
% 
% mesh.phase.colpoints = 4*ones(1,10)
% mesh.phase(3).colpoints = 5*ones(1,20)
% mesh.phase(4).colpoints = 4*ones(1,10)
%--------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%

setup.name                           = 'SPARTAN-Combined';
setup.functions.continuous           = @CombinedContinuous;
setup.functions.endpoint             = @CombinedEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';

setup.nlp.ipoptoptions.maxiterations = 350;

setup.derivatives.supplier           = 'sparseFD';

setup.derivatives.derivativelevel    = 'first';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-guessUpdate';
setup.derivatives.dependencies      = 'full';



%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%

%% Parallel Loop

num_it = 3; % Define number of iterations

% Create variable setup structure
for i = 1:num_it
setup_par(i) = setup;
% setup_par(i).nlp.ipoptoptions.maxiterations = 500 + 10*i;
setup_par(i).guess.phase.state(:,1)   = [24000;33000 + 1000*i]; % vary altitude guess
end

parfor i = 1:num_it
    
output_temp = gpops2(setup_par(i)); % Run GPOPS-2. Use setup for each parallel iteration.

% Run forward simulation for comparison between runs
% Extract states
% alt22 = output_temp.result.solution.phase(3).state(:,1).';
% lon22 = output_temp.result.solution.phase(3).state(:,2).';
% lat22 = output_temp.result.solution.phase(3).state(:,3).';
% v22 = output_temp.result.solution.phase(3).state(:,4).'; 
% gamma22 = output_temp.result.solution.phase(3).state(:,5).'; 
% zeta22 = output_temp.result.solution.phase(3).state(:,6).';
% alpha22 = output_temp.result.solution.phase(3).state(:,7).';
% eta22 = output_temp.result.solution.phase(3).state(:,8).';
% mFuel22 = output_temp.result.solution.phase(3).state(:,9).'; 
% time22 = output_temp.result.solution.phase(3).time.';
% throttle22 = output_temp.result.solution.phase(3).state(:,10).';
% 
% % forward simulation
% forward0 = [alt22(1),gamma22(1),v22(1),zeta22(1),lat22(1),lon22(1), mFuel22(1)];
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1):time22(end),forward0);

% error(i) = (f_y(end,6) + lon2(end))^2 + (f_y(end,5) + lat2(end))^2;
% error(i) = abs(mFuel22(end) - f_y(end,7));
error(i) = output_temp.result.objective;

output_store{i} = output_temp;

end

[min_error,index] = min(error); % Calculate the result which minimises the chosen error function

% [max_pl,index] = max(PayloadMass);% Calculate the result which maximises payload mass the chosen error function
output = output_store{index};


%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================
% Assign the primal variables
alt21 = output.result.solution.phase.state(:,1);

lon21 = output.result.solution.phase.state(:,2);

lat21 = output.result.solution.phase.state(:,3);
v21 = output.result.solution.phase.state(:,4); 
gamma21 = output.result.solution.phase.state(:,5); 
zeta21 = output.result.solution.phase.state(:,6);
alpha21 = output.result.solution.phase.state(:,7);
eta21 = output.result.solution.phase.state(:,8);
mFuel21 = output.result.solution.phase.state(:,9); 

aoadot21  = output.result.solution.phase.control(:,1); 
etadot21  = output.result.solution.phase.control(:,2); 


time21 = output.result.solution.phase.time;
[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho21,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,1, mFuel21,mFuel21(1),mFuel21(end), 1, 0);

%% plot Ascent
figure(211)
fig = gcf;
set(fig,'Position',[200 0 850 1200])

subplot(7,2,1)
hold on
plot(time21, alt21/1000,'Color','k')
title('Trajectory (km)')

subplot(7,2,3)
hold on
plot(time21, v21,'Color','k')
title('Velocity (m/s)')

subplot(7,2,4)
plot(time21, M21,'Color','k')
title('Mach no')

subplot(7,2,5)
plot(time21, q21/1000,'Color','k')
title('Dynamic Pressure (kpa)')

subplot(7,2,6)
hold on
plot(time21, rad2deg(gamma21),'Color','k')
title('Trajectory Angle (Deg)')


subplot(7,2,7)

plot(time21, rad2deg(alpha21),'Color','k')
title('Angle of Attack (deg)')

subplot(7,2,8)
hold on

plot(time21, rad2deg(eta21),'Color','k')
title('Bank Angle (deg)')

subplot(7,2,9)
plot(time21, flapdeflection21,'Color','k')
title('Flap Deflection (deg)')


% Isp1 = T1./Fueldt1./9.81;
IspNet1 = (T21-Fd21)./Fueldt21./9.81;

subplot(7,2,10)
plot(time21, T21/1000,'Color','k')
title('Thrust (kN)')

subplot(7,2,11)
plot(time21, IspNet1,'Color','k')
title('Net Isp (s)')
xlabel('Time (s)');

subplot(7,2,12)
plot(time21, mFuel21,'Color','k')
title('Fuel Mass (kg)')
xlabel('Time (s)');

subplot(7,2,13)
plot(time21, rad2deg(zeta21),'Color','k')
title('Heading Angle (deg)')
xlabel('Time (s)');


SecondStageStates = [[time21 time22]' [alt21 alt22]' [lon21 lon22]' [lat21 lat22]' [v21 v22]' [gamma21 gamma22]' [zeta21 zeta22]' [alpha21 alpha22]' [eta21 eta22]' [mFuel21 mFuel22]'];
dlmwrite('SecondStageStates',['time (s) ' 'altitude (m) ' 'longitude (rad) ' 'latitude (rad) ' 'velocity (m/s) ' 'trajectory angle (rad) ' 'heading angle (rad) ' 'angle of attack (rad) ' 'bank angle (rad) ' 'fuel mass (kg) '],'');
dlmwrite('SecondStageStates',SecondStageStates,'-append','delimiter',' ');
copyfile('SecondStageStates',sprintf('../ArchivedResults/%s/SecondStage_%s',strcat(Timestamp,'mode',num2str(mode)),Timestamp));

%% Plot Return
figure(221)
fig = gcf;
set(fig,'Position',[200 0 850 1200])

subplot(7,2,1)
xlim([time22(1) time22(end)]);
hold on
plot(time22, alt22/1000,'Color','k')
title('Trajectory (km)')

dim = [.55 .7 .2 .2];
annotation('textbox',dim,'string',{['Fuel Used: ' num2str(mFuel22(1)) ' kg']},'FitBoxToText','on');  

subplot(7,2,2)
xlim([time22(1) time22(end)]);
hold on
plot(time22, v22,'Color','k')
title('Velocity (m/s)')

subplot(7,2,3)
hold on
xlim([time22(1) time22(end)]);
plot(time22, M22,'Color','k')
title('Mach no')

subplot(7,2,4)
hold on
xlim([time22(1) time22(end)]);
plot(time22, q22/1000,'Color','k')
title('Dynamic Pressure (kpa)')

subplot(7,2,5)
xlim([time22(1) time22(end)]);
hold on
plot(time22, rad2deg(gamma22),'Color','k')
title('Trajectory Angle (Deg)')

subplot(7,2,6)
hold on
xlim([time22(1) time22(end)]);
plot(time22, rad2deg(alpha22),'Color','k')
title('Angle of Attack (deg)')

subplot(7,2,7)
xlim([time22(1) time22(end)]);
ylim([rad2deg(min(eta22))-1 rad2deg(max(eta22))+1])
hold on
plot(time22, rad2deg(eta22),'Color','k')
title('Bank Angle (deg)')

subplot(7,2,8)
hold on
xlim([time22(1) time22(end)]);
plot(time22, flapdeflection22,'Color','k')
title('Flap Deflection (deg)')

subplot(7,2,9)
hold on
xlim([time22(1) time22(end)]);
plot(time22, T22/1000,'Color','k')
title('Thrust (kN)')

subplot(7,2,10)
hold on
xlim([time22(1) time22(end)]);
plot(time22, Isp22,'Color','k')
title('Potential Isp (s)')

subplot(7,2,11)
hold on
xlim([time22(1) time22(end)]);
plot(time22, mFuel22,'Color','k')
title('Fuel Mass (kg)')
xlabel('Time (s)');

subplot(7,2,12)
hold on
xlim([time22(1) time22(end)]);
plot(time22, throttle22,'Color','k')
title('Throttle')
xlabel('Time (s)');

subplot(7,2,13)
hold on
xlim([time22(1) time22(end)]);
plot(time22, rad2deg(zeta22),'Color','k')
title('Heading Angle (deg)')
xlabel('Time (s)');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt21(1),gamma21(1),v21(1),zeta21(1),lat21(1),lon21(1), mFuel21(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time21,alpha21,f_t),ControlInterp(time21,eta21,f_t),1,mFuel21(1),mFuel21(end)),time21(1:end),forward0);

figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time21,alt21);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time21,gamma21);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time21,v21);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time21,zeta21);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time21,mFuel21);



% Return Forward
forward0 = [alt22(1),gamma22(1),v22(1),zeta22(1),lat22(1),lon22(1), mFuel22(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1):time22(end),forward0);

figure(213)
subplot(7,1,1)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time22,alt22);

% gamma  = output.result.solution.phase.state(:,5);

subplot(7,1,2)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time22,gamma22);

% latitude  = output.result.solution.phase.state(:,3);
subplot(7,1,3:5)
hold on
plot(f_y(:,6),f_y(:,5));
plot(lon22,lat22);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time22,mFuel22);

%% Check KKT and pontryagins minimum
% Check that the hamiltonian = 0 (for free end time)
% Necessary condition
input_test = output.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

lambda1 = output.result.solution.phase.costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test.dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output.result.solution.phase(3).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(3).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

figure(2410)
hold on
plot(time21(1:end-1),H1)
plot(time22(1:end-1),H2)
ylabel('Hamiltonian')
xlabel('Time (s)')
legend('Ascent','Return')

%% Check State Feasibility
% Check calculated derivatives with the numerical derivative of each
% porimal, scaled by that primal
figure(2420)
hold on
for i = 1:length(output.result.solution.phase.state(1,:))
plot(time21,([diff(output.result.solution.phase.state(:,i))./diff(output.result.solution.phase.time); 0] - phaseout_test.dynamics(:,i))./output.result.solution.phase.state(:,i),'--');
end
for i = 1:length(output.result.solution.phase(3).state(1,:))
    if i<= 7 % Plot different line styles when no. of colours exceeded
    plot(time22,([diff(output.result.solution.phase(3).state(:,i))./diff(output.result.solution.phase(3).time); 0] - phaseout_test(3).dynamics(:,i))./output.result.solution.phase(3).state(:,i));
    else
    plot(time22,([diff(output.result.solution.phase(3).state(:,i))./diff(output.result.solution.phase(3).time); 0] - phaseout_test(3).dynamics(:,i))./output.result.solution.phase(3).state(:,i),':');
    end
end
title('Derivative Feasibility Check')
xlabel('Time (s)')
ylabel('Derivative Error')
ylim([-1,1])
legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')

%% plot engine interpolation visualiser
T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt21); 
T_in1 = auxdata.interp.tempgridded(M21,rad2deg(alpha21)).*T0;
M_in1 = auxdata.interp.M1gridded(M21, rad2deg(alpha21));

plotM = [min(M_englist):0.01:9];
plotT = [min(T_englist):1:550];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = auxdata.interp.eqGridded(gridM,gridT);
interpIsp = auxdata.interp.IspGridded(gridM,gridT);

figure(2100)
hold on
contourf(gridM,gridT,interpeq,100,'LineWidth',0.0);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'k');
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'Equivalence Ratio';
caxis([.4 1])
plot(M_in1,T_in1,'r');

error_Isp = auxdata.interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);

figure(2110)
hold on
contourf(gridM,gridT,interpIsp,100,'LineWidth',0);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'ISP';
plot(M_in1,T_in1,'r');

figure(2120)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.M1gridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'M1';

figure(2130)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.tempgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'T1/T0';

figure(2140)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.presgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'P1/P0';

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);


%% ThirdStage

alt3  = output.result.solution.phase(4).state(:,1);
v3    = output.result.solution.phase(4).state(:,2);
gamma3  = output.result.solution.phase(4).state(:,3);
m3    = output.result.solution.phase(4).state(:,4);
aoa3    = output.result.solution.phase(4).state(:,5);
phi3    = output.result.solution.phase(4).state(:,6);
zeta3    = output.result.solution.phase(4).state(:,7);
aoadot3       = output.result.solution.phase(4).control(:,1);

forward0 = [alt3(1),v3(1),gamma3(1),m3(1),phi3(1),zeta3(1)];

time3 = output.result.solution.phase(4).time;

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel3_forward(f_t, f_y,auxdata,ControlInterp(time3,aoa3,f_t),ControlInterp(time3,aoadot3,f_t)),time3(1:end),forward0);

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDyn(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,phi3,zeta3);


xi3(1) = lon21(end);
for i = 2:length(time3)
    xi3(i) = xi3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end

[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,phiexo,incexo,Texo,CLexo,Lexo,incdiffexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), phi3(end),xi3(end), zeta3(end), m3(end), auxdata);




figure(312)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time3,alt3);

figure(313)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time3,v3);



figure(311)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
    
    subplot(4,2,1)
    hold on
    title('Altitude (km');
    plot([time3; timeexo.'+time3(end)], [alt3; altexo.']/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,2)
    hold on
    title('Dynamic Pressure (kPa');
    plot([time3; timeexo.'+time3(end)],[q3;qexo.';qexo(end)]/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,3)
    hold on
    title('Angle of Attack (deg)');
    plot([time3; timeexo.'+time3(end)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,4)
    hold on
    title('Velocity (m/s)');
    plot([time3; timeexo.'+time3(end)],[v3;v3exo.'],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,5)
    hold on
    title('Mass (kg');
    plot([time3; timeexo.'+time3(end)],[ m3;mexo.';mexo(end)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,6)
    hold on
    title('Thrust Vector Angle (deg)');
    plot([time3; timeexo.'+time3(end)],[rad2deg(Vec_angle3);0*ones(length(timeexo),1)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,7)
    hold on
    title('Trajectory Angle (deg)');
    plot([time3; timeexo.'+time3(end)], [rad2deg(gamma3);rad2deg(gammaexo).'],'Color','k')

    xlabel('Time (s)');
    xlim([time3(1) timeexo(end)+time3(end)])

    
    % Write data to file
    dlmwrite('ThirdStageData',['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg) ' 'dynamic pressure (Pa)' 'trajectory angle (rad) ' 'Lift (N)' 'Drag (N)' 'heading angle (rad) ' 'latitude (rad) ' 'angle of attack (rad) '],'');
    dlmwrite('ThirdStageData',[[time3; time3(end)+timeexo'], [alt3; altexo'], [v3; v3exo'], [m3; mexo'; mexo(end)],[q3; qexo'; qexo(end)] ,[gamma3; gammaexo'],[L3; Lexo'; Lexo(end)],[D3; Dexo'; Dexo(end)] ,[zeta3; zetaexo'], [phi3; phiexo'], [aoa3; zeros(length(timeexo),1)]],'-append','delimiter',' ')
copyfile('ThirdStageData',sprintf('../ArchivedResults/%s/ThirdStage_%s',strcat(Timestamp,'mode',num2str(mode)),Timestamp));


%% SAVE FIGS

saveas(figure(311),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'ThirdStage.fig']);
print(figure(311),'ThirdStage','-dpng');
movefile('ThirdStage.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(211),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'SecondStage.fig']);
saveas(figure(221),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'Return.fig']);    
saveas(figure(2410),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'Hamiltonian.fig']);
saveas(figure(2420),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'Validation.fig']);
saveas(figure(212),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'Forward1.fig']);
saveas(figure(213),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'Forward2.fig']);
% saveas(figure(2100),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'eq.fig']);
% saveas(figure(2110),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ISP.fig']);
saveas(figure(2301),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'GroundTrack.fig']);

% 
%% First Stage =========================================================

t1 = output.result.solution.phase(1).time.';
% 
alt1 = output.result.solution.phase(1).state(:,1).';
v1 = output.result.solution.phase(1).state(:,2).';
m1 = output.result.solution.phase(1).state(:,3).';
gamma1 = output.result.solution.phase(1).state(:,4).';
alpha1 = output.result.solution.phase(1).state(:,5).';
zeta1 = output.result.solution.phase(1).state(:,6).';
phi1 = output.result.solution.phase(1).state(:,8).';
xi1 = output.result.solution.phase(1).state(:,9).';
% 

% 
FirstStageSMF = (mRocket - mFuel)/(m1(1) - mSpartan);
% 

FirstStageStates = [t1' alt1' v1' m1' gamma1' alpha1' zeta1' phi1' xi1'];

dlmwrite('FirstStageStates',['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg)' 'trajectory angle (rad) ' 'angle of attack (rad) ' 'heading angle (rad) ' 'latitude (rad)'],'');
dlmwrite('FirstStageStates',FirstStageStates,'-append','delimiter',' ');
copyfile('FirstStageStates',sprintf('../ArchivedResults/%s/FirstStage_%s',strcat(Timestamp,'mode',num2str(mode)),Timestamp));


% Iterative Prepitch Determination ========================================
%This back determines the mass and launch altitude necessary to get to
%100m, 30m/s at the PS method determined fuel mass


interp = auxdata.interp;
Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

% ntoe that launch altitude does vary, but it should only be slightly
controls = fminunc(@(controls) prepitch(controls,m1(1),interp,Throttle,Vehicle,Atmosphere),[10,6]);


h_launch = controls(1)
t_prepitch = controls(2)
Isp1 = Vehicle.Isp.SL;
T1 = Vehicle.T.SL;
dm1 = -T1./Isp1./9.81;
m0_prepitch = m1(1) - dm1*t_prepitch;



%% Forward Integrator
 phase = 'postpitch';
tspan = t1; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [h0 v0 m1(1) deg2rad(89.9) phi1(1) zeta1(1)];

[t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t_f,t1,zeta1),ControlFunction(t_f,t1,alpha1),phase,interp,Throttle,Vehicle,Atmosphere), tspan, postpitch0_f);

figure(103)
hold on
plot(postpitch_f(:,1));
plot(alt1);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = h_launch;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan2 = [0 t_prepitch]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0, 0, 0];

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,interp,Throttle,Vehicle,Atmosphere), tspan2, y0);  


figure(111);
hold on
title('First Stage Trajectory');
    fig = gcf;
set(fig,'Position',[200 0 850 600])
subplot(4,2,1)
hold on
title('Trajectory Angle (deg)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [rad2deg(y(:,4).') rad2deg(gamma1)],'color','k');
subplot(4,2,2)
hold on
title('Velocity (m/s)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [y(:,2).' v1],'color','k');
subplot(4,2,3)
hold on
title('Altitude (km)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [y(:,1).'/1000 alt1/1000],'color','k');
subplot(4,2,4)
hold on
title('Angle of Attack (deg)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha1)],'color','k');
subplot(4,2,5)
hold on
title('Mass (kg)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [y(:,3).' m1],'color','k');
subplot(4,2,6)
hold on
title('Heading Angle (deg)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [rad2deg(y(:,6).') rad2deg(zeta1)],'color','k');
subplot(4,2,7)
hold on
title('Latitude (deg)');
xlim([0 t1(end)+t_prepitch(end)]);
plot([t_prepitch.' t1+t_prepitch(end)], [rad2deg(phi1(1)+y(:,8).') rad2deg(phi1)],'color','k');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')
xlabel('Time (s)')
xlim([0,t1(end)+t_prepitch(end)]);

saveas(figure(111),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,'FirstStage.fig']);



%% Create Easy Latex Inputs

dlmwrite('LatexInputs.txt',strcat('\newcommand{\PayloadToOrbitMode', num2str(mode) ,'}{ ', num2str(round(ThirdStagePayloadMass,1),'%.1f') , '}'), 'delimiter','','newline', 'pc')
dlmwrite('LatexInputs.txt',strcat('\newcommand{\12SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(1)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc')

dlmwrite('LatexInputs.txt',strcat('\newcommand{\FirstStageSMFMode', num2str(mode) ,'}{ ', num2str(round(FirstStageSMF,3),'%.3f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite('LatexInputs.txt',strcat('\newcommand{\23SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(end)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite('LatexInputs.txt',strcat('\newcommand{\23SeparationvMode', num2str(mode) ,'}{ ', num2str(round(v21(end),0)) , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite('LatexInputs.txt',strcat('\newcommand{\23SeparationqMode', num2str(mode) ,'}{ ', num2str(round(q21(end)/1000,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite('LatexInputs.txt',strcat('\newcommand{\23SeparationLDMode', num2str(mode) ,'}{ ', num2str(round(L21(end)/Fd21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite('LatexInputs.txt',strcat('\newcommand{\2FlightTimeMode', num2str(mode) ,'}{ ', num2str(round(time21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

qlt20 = find(q3<20000);
dlmwrite('LatexInputs.txt',strcat('\newcommand{\3qOver20Mode', num2str(mode) ,'}{ ', num2str(round(time3(qlt20(1))-time3(1),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');



%% Bound Check
% Peform check to see if any of the states are hitting their bounds. This
% is an error if the bound is not intended to constrain the state. Fuel
% mass and throttle are not checked, as these will always hit bounds. 

for i = 1: length(output.result.solution.phase.state(1,:))
    if any(output.result.solution.phase(1).state(:,i) == bounds.phase(1).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting lower bound'))
    end
    
    if any(output.result.solution.phase(1).state(:,i) == bounds.phase(1).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting upper bound'))
    end
end

for i = 1: length(output.result.solution.phase.state(1,:))-1
    if any(output.result.solution.phase.state(:,i) == bounds.phase.state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting lower bound'))
    end
    
    if any(output.result.solution.phase.state(:,i) == bounds.phase.state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting upper bound'))
    end
end

for i = 1: length(output.result.solution.phase(3).state(1,:))-2
    if any(output.result.solution.phase(3).state(:,i) == bounds.phase(3).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting lower bound'))
    end
    
    if any(output.result.solution.phase(3).state(:,i) == bounds.phase(3).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting upper bound'))
    end
end

% Angle of attack is not checked on third stage, because angle of attack is hard constrained and should be checked manually. 
for i = [1:3 6: length(output.result.solution.phase(4).state(1,:))]
    if any(output.result.solution.phase(4).state(:,i) == bounds.phase(4).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting lower bound'))
    end
    
    if any(output.result.solution.phase(4).state(:,i) == bounds.phase(4).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting upper bound'))
    end
end

%% =========================================================================
% Troubleshooting Procedure
% =========================================================================
% 1: Check that you have posed your problem correctly ie. it is physically
% feasible and the bounds allow for a solution.
% 2: Check if there is any extrapolations causing bad dynamics. 
% 3: Check guess, is it reasonable?.
% 4: Increase number of iterations and collocation points.
% 5: Check for large nonlinearities, eg. atmospheric properties suddenly going to zero or thrust cutting off. These need to be eliminated or separated into phases.
% 6: Check for NaN values (check derivatives in Dynamics file while
% running).






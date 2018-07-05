%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc


auxdata.delta = deg2rad(0) % thrust vector angle test
auxdata.dragmod = 1. %drag increase test


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\..\thirdStage')
addpath('..\EngineData')
addpath('..\')

Timestamp = datestr(now,30)
mkdir('../ArchivedResults', sprintf(Timestamp))
copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',Timestamp))

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
lonMin = -pi;         lonMax = -lonMin;
latMin = -70*pi/180;  latMax = -latMin;
lat0 = -0.264;
lon0 = deg2rad(145);
aoaMin = 0;  aoaMax = 10*pi/180;
bankMin1 = -1*pi/180; bankMax1 =   90*pi/180;

% Primal Bounds
% bounds.phase(1).state.lower = [Stage2.Bounds.Alt(1), lonMin, latMin, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin, bankMin1, Stage2.Bounds.mFuel(1)];
% bounds.phase(1).state.upper = [Stage2.Bounds.Alt(2), lonMax, latMax, Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax, bankMax1, Stage2.Bounds.mFuel(2)];
bounds.phase(1).state.lower = [Stage2.Bounds.Alt(1), 2, -0.5, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin, bankMin1, Stage2.Bounds.mFuel(1)];
bounds.phase(1).state.upper = [Stage2.Bounds.Alt(2), 3, 0.5, Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax, bankMax1, Stage2.Bounds.mFuel(2)];

% Initial States
bounds.phase(1).initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, Stage2.Initial.v, Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin, bankMin1, Stage2.Initial.mFuel] ;
bounds.phase(1).initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, Stage2.Initial.v, Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax, bankMax1, Stage2.Initial.mFuel];

% End States

% bounds.phase(1).finalstate.lower = [Stage2.Bounds.Alt(1), lonMin, latMin, Stage2.Bounds.v(1), Stage2.End.gammaOpt(1), Stage2.End.Zeta, aoaMin, 0, Stage2.End.mFuel];
% bounds.phase(1).finalstate.upper = [Stage2.Bounds.Alt(2), lonMax, latMax, Stage2.Bounds.v(2), Stage2.End.gammaOpt(2), Stage2.End.Zeta, aoaMax, 0, Stage2.Initial.mFuel];
bounds.phase(1).finalstate.lower = [34000, 2, -0.5, 2500, 0, Stage2.Bounds.zeta(1), aoaMin, 0, Stage2.End.mFuel];
bounds.phase(1).finalstate.upper = [45000, 3, 0.5, Stage2.Bounds.v(2), deg2rad(20), Stage2.Bounds.zeta(2), aoaMax, 0, Stage2.Initial.mFuel];

% Control Bounds
% bounds.phase(1).control.lower = [deg2rad(-.1), deg2rad(-1)];
% bounds.phase(1).control.upper = [deg2rad(.1), deg2rad(1)];
bounds.phase(1).control.lower = [deg2rad(-.1), deg2rad(-1)];
bounds.phase(1).control.upper = [deg2rad(.1), deg2rad(1)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = Stage2.Bounds.time(1);
bounds.phase(1).finaltime.upper = Stage2.Bounds.time(2);

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.
if mode == 1 || mode == 14 || mode == 15
    bounds.phase(1).path.lower = [0];
    bounds.phase(1).path.upper = [50000];
elseif mode == 12
    bounds.phase(1).path.lower = [0];
    bounds.phase(1).path.upper = [55000];
elseif mode == 13
    bounds.phase(1).path.lower = [0];
    bounds.phase(1).path.upper = [45000];
elseif mode ==3 || mode == 32
        bounds.phase(1).path.lower = [0];
    bounds.phase(1).path.upper = [50010];
end


%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
% guess.phase(1).state(:,1)   = [24000;35000];
guess.phase(1).state(:,2)   = [2.53;2.5368];
guess.phase(1).state(:,3)   = [-0.269;-0.10];
guess.phase(1).state(:,4)   = Stage2.Guess.v.';
guess.phase(1).state(:,5)   = Stage2.Guess.gamma.';
guess.phase(1).state(:,6)   = Stage2.Guess.zeta.';
guess.phase(1).state(:,7)   = [2*pi/180; 5*pi/180];
guess.phase(1).state(:,8)   = [deg2rad(10);deg2rad(10)];
guess.phase(1).state(:,9) 	= [Stage2.Initial.mFuel, 100];

guess.phase(1).control      = [[0;0],[0;0]];
guess.phase(1).time          = [0;650];

% Tie stages together
bounds.eventgroup(1).lower = [zeros(1,10)];
bounds.eventgroup(1).upper = [zeros(1,10)]; 

%% Flyback
tfMin = 0;            tfMax = 5000;
altMin = 10;  altMax = 70000;
speedMin = 10;        speedMax = 5000;
fpaMin = -80*pi/180;  fpaMax =  80*pi/180;
aziMin = 60*pi/180; aziMax =  500*pi/180;
mFuelMin = 0; mFuelMax = 500;
bankMin2 = -10*pi/180; bankMax2 =   90*pi/180;

% lonf = deg2rad(145);
% latf   = -0.269;

throttleMin = 0; throttleMax = 1;

bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = 3000;
bounds.phase(2).finaltime.lower = 400;
bounds.phase(2).finaltime.upper = 4000;

% Initial Bounds
% bounds.phase(2).initialstate.lower = [altMin, lonMin, latMin, speedMin, fpaMin, aziMin, aoaMin, bankMin2, mFuelMin, throttleMin];
% bounds.phase(2).initialstate.upper = [altMax, lonMax, latMax, speedMax, fpaMax, aziMax, aoaMax, bankMax2, mFuelMax, throttleMax];

bounds.phase(2).initialstate.lower = [altMin, lonMin, latMin, speedMin, fpaMin, aziMin, aoaMin, 0, mFuelMin, throttleMin];
bounds.phase(2).initialstate.upper = [altMax, lonMax, latMax, speedMax, fpaMax, aziMax, aoaMax, 0, mFuelMax, throttleMax];

% State Bounds
bounds.phase(2).state.lower = [altMin, lonMin, latMin, speedMin, fpaMin, aziMin, aoaMin, bankMin2, mFuelMin, throttleMin];
bounds.phase(2).state.upper = [altMax, lonMax, latMax, speedMax, fpaMax, aziMax, aoaMax, bankMax2, mFuelMax, throttleMax];

% End State Bounds
% bounds.phase(2).finalstate.lower = [altMin, lonf-0.002, latf-0.002, speedMin, deg2rad(-10), aziMin, aoaMin, bankMin2, Stage2.End.mFuel, throttleMin];
% bounds.phase(2).finalstate.upper = [200, lonf+0.002, latf+0.002, speedMax, deg2rad(30), aziMax, aoaMax, bankMax2, Stage2.End.mFuel, throttleMax];
bounds.phase(2).finalstate.lower = [altMin, lonMin-0.001, latMin-0.001, speedMin, deg2rad(-20), aziMin, aoaMin, bankMin2, Stage2.End.mFuel, throttleMin];
% bounds.phase(2).finalstate.upper = [500, lonMax+0.001, latMax+0.001, speedMax, deg2rad(30), aziMax, aoaMax, bankMax2, Stage2.End.mFuel, throttleMax];
bounds.phase(2).finalstate.upper = [altMax, lonMax+0.001, latMax+0.001, speedMax, deg2rad(30), aziMax, aoaMax, bankMax2, Stage2.End.mFuel, throttleMax];

% Control Bounds
bounds.phase(2).control.lower = [deg2rad(-.2), deg2rad(-5), -1];
bounds.phase(2).control.upper = [deg2rad(.2), deg2rad(5), 1];

% Path Bounds
bounds.phase(2).path.lower = 0;
bounds.phase(2).path.upper = 50000;

bounds.eventgroup(2).lower = [-0.001 -0.001]; % Constrain final latitude and longitude, with variable first stage approximation
bounds.eventgroup(2).upper = [0.001 0.001]; 

% Guess
tGuess              = [440; 1500];
altGuess            = [35000; 100];
lonGuess            = [lon0; lon0-.1*pi/180];
latGuess            = [-0.11;-0.10-0*pi/180];
speedGuess          = [3000; 10];
fpaGuess            = [0; 0];
aziGuess            = [deg2rad(97); deg2rad(270)];
aoaGuess            = [8*pi/180; 5*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];
% mFuelGuess          = [mFuelMax; mFuelMin];
mFuelGuess          = [100; mFuelMin];
guess.phase(2).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[0.;0.]];
guess.phase(2).control = [[0;0],[0;0],[0;0]];
% guess.phase.control = [aoaGuess];
guess.phase(2).time    = tGuess;


%% Third Stage

tfMin3 = 7;            tfMax3 = 200;
altMin3 = 30000;  altMax3 = 84000;
phiMin3 = -0.2;         phiMax3 = 0.1;
vMin3 = 10;        vMax3 = 8000;
gammaMin3=deg2rad(-5);  gammaMax3 =  deg2rad(30);
zetaMin3 = deg2rad(90); zetaMax3 =  deg2rad(110);
aoaMin3 = 0;  aoaMax3= deg2rad(20);

aoadotMin3 = -deg2rad(1);
aoadotMax3 = deg2rad(1);

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase(3).initialtime.lower = 0;
bounds.phase(3).initialtime.upper = 10000;
bounds.phase(3).finaltime.lower = 1;
bounds.phase(3).finaltime.upper = 10000;


% bounds.phase.initialstate.lower = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMin, phi0, zeta0];
% bounds.phase.initialstate.upper = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMax, phi0, zeta0];
% bounds.phase.initialstate.lower = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMin, phi0, zetaMin];
% bounds.phase.initialstate.upper = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMax, phi0, zetaMax];

bounds.phase(3).initialstate.lower = [altMin3,vMin3, deg2rad(1),  auxdata.Stage3.mTot, aoaMin3, phiMin3, zetaMin3];
bounds.phase(3).initialstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];


% bounds.phase(3).initialstate.lower = [altMin3,vMin3, deg2rad(1),  2000, aoaMin3, phiMin3, zetaMin3];
% bounds.phase(3).initialstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(3).state.lower = [altMin3,vMin3, gammaMin3, 0, aoaMin3, phiMin3, zetaMin3];
bounds.phase(3).state.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(3).finalstate.lower = [altMin3, vMin3, 0, 0, 0, phiMin3, zetaMin3];
bounds.phase(3).finalstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, 0, phiMax3, zetaMax3];

bounds.phase(3).control.lower = [aoadotMin3];
bounds.phase(3).control.upper = [aoadotMax3];

bounds.phase(3).path.lower = [-deg2rad(8), -inf];
bounds.phase(3).path.upper = [deg2rad(8), 0];


bounds.eventgroup(3).lower = [zeros(1,7) 90000 0];
bounds.eventgroup(3).upper = [zeros(1,6) 1000 566000 0];

tGuess              = [0; 150];
altGuess            = [35000; 60000];
vGuess          = [2700; 6000];
gammaGuess            = [0; deg2rad(10)];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(20); deg2rad(20)];
phiGuess = [-0.11;-0.11];
zetaGuess = [deg2rad(97);deg2rad(97)];
guess.phase(3).state   = [altGuess, vGuess, gammaGuess, mGuess, aoaGuess, phiGuess, zetaGuess];
guess.phase(3).control = [0;0];
guess.phase(3).time    = tGuess;

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre'; % Default method does not perform h adaptions sometimes, and can be unstable. 
%  mesh.method       = 'hp-DarbyRao';
mesh.maxiterations = 2;
mesh.colpointsmin = 3;
mesh.colpointsmax = 50;
mesh.tolerance    = 1e-5;

% mesh.phase(1).fraction = 0.1*ones(1,10)
% 
% mesh.phase(2).fraction = 0.05*ones(1,20)
% 
% mesh.phase(3).fraction = 0.1*ones(1,10)
% 
% mesh.phase(1).colpoints = 4*ones(1,10)
% mesh.phase(2).colpoints = 5*ones(1,20)
% mesh.phase(3).colpoints = 4*ones(1,10)
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

setup.nlp.ipoptoptions.maxiterations = 500;

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
setup_par(i).guess.phase(1).state(:,1)   = [24000;33000 + 1000*i]; % vary altitude guess
end

parfor i = 1:num_it
    
output_temp = gpops2(setup_par(i)); % Run GPOPS-2. Use setup for each parallel iteration.

% Run forward simulation for comparison between runs
% Extract states
alt2 = output_temp.result.solution.phase(2).state(:,1).';
lon2 = output_temp.result.solution.phase(2).state(:,2).';
lat2 = output_temp.result.solution.phase(2).state(:,3).';
v2 = output_temp.result.solution.phase(2).state(:,4).'; 
gamma2 = output_temp.result.solution.phase(2).state(:,5).'; 
zeta2 = output_temp.result.solution.phase(2).state(:,6).';
Alpha2 = output_temp.result.solution.phase(2).state(:,7).';
eta2 = output_temp.result.solution.phase(2).state(:,8).';
mFuel2 = output_temp.result.solution.phase(2).state(:,9).'; 
time2 = output_temp.result.solution.phase(2).time.';
throttle2 = output_temp.result.solution.phase(2).state(:,10).';

% forward simulation
forward0 = [alt2(1),gamma2(1),v2(1),zeta2(1),lat2(1),lon2(1), mFuel2(1)];
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time2,Alpha2,f_t),ControlInterp(time2,eta2,f_t),ThrottleInterp(time2,throttle2,f_t)),time2(1):time2(end),forward0);

% error(i) = (f_y(end,6) + lon2(end))^2 + (f_y(end,5) + lat2(end))^2;
error(i) = abs(mFuel2(end) - f_y(end,7));
output_store{i} = output_temp;

end

[min_error,index] = min(error); % Calculate the result which minimises the chosen error function

output = output_store{index};


%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================
% Assign the primal variables
alt = output.result.solution.phase(1).state(:,1).';
alt2 = output.result.solution.phase(2).state(:,1).';
lon = output.result.solution.phase(1).state(:,2).';
lon2 = output.result.solution.phase(2).state(:,2).';
lat = output.result.solution.phase(1).state(:,3).';
lat2 = output.result.solution.phase(2).state(:,3).';
v = output.result.solution.phase(1).state(:,4).'; 
v2 = output.result.solution.phase(2).state(:,4).'; 
gamma = output.result.solution.phase(1).state(:,5).'; 
gamma2 = output.result.solution.phase(2).state(:,5).'; 
zeta = output.result.solution.phase(1).state(:,6).';
zeta2 = output.result.solution.phase(2).state(:,6).';
Alpha = output.result.solution.phase(1).state(:,7).';
Alpha2 = output.result.solution.phase(2).state(:,7).';
eta = output.result.solution.phase(1).state(:,8).';
eta2 = output.result.solution.phase(2).state(:,8).';
mFuel = output.result.solution.phase(1).state(:,9).'; 
mFuel2 = output.result.solution.phase(2).state(:,9).'; 

throttle2 = output.result.solution.phase(2).state(:,10).';

aoadot1  = output.result.solution.phase(1).control(:,1).'; 
etadot1  = output.result.solution.phase(1).control(:,2).'; 

aoadot2  = output.result.solution.phase(2).control(:,1).'; 
etadot2  = output.result.solution.phase(2).control(:,2).'; 

time = output.result.solution.phase(1).time.';
time2 = output.result.solution.phase(2).time.';

figure(201)
subplot(9,1,1)
hold on
plot(time,alt)
plot(time2,alt2)
subplot(9,1,2)
hold on
plot(time,v)
plot(time2,v2)
subplot(9,1,3)
hold on
plot(time,lon)
plot(time2,lon2)
subplot(9,1,4)
hold on
plot(time,lat)
plot(time2,lat2)
subplot(9,1,5)
hold on
plot(time,v)
plot(time2,v2)
subplot(9,1,6)
hold on
plot(time,gamma)
plot(time2,gamma2)
subplot(9,1,7)
hold on
plot(time,ones(1,length(time)))
plot(time2,throttle2)


figure(230)
hold on
plot3(lon,lat,alt)
plot3(lon2,lat2,alt2)

  
    figure(2301)
hold on

axesm('pcarree','Origin',[0 rad2deg(lon0) 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
% plotm(rad2deg(lat),rad2deg(lon+lon0))
plotm(rad2deg(lat),rad2deg(lon),'b')
plotm(rad2deg(lat2),rad2deg(lon2),'r')
    
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)

% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point

global phi

% cd('../ThirdStage')
% [ThirdStagePayloadMass,ThirdStageControls,ThirdStageZeta,ThirdStagePhi,ThirdStageAlt,ThirdStagev,ThirdStaget,ThirdStageAlpha,ThirdStagem,ThirdStagegamma,ThirdStageq] = ThirdStageOptm(V(end),gamma(end),v(end), phi(end),zeta(end), 1);
% ThirdStagePayloadMass
% cd('../SecondStage')
ThirdStagePayloadMass = 0;
% 
% 
% if auxdata.PayloadGrid(alt(end)+10,gamma(end),v(end)) - auxdata.PayloadGrid(alt(end),gamma(end),v(end)) < 0
%     disp('Check Third Stage Payload Matrix, Found Maxima')
% end
% 
% if gamma(end) > max(ThirdStageData(:,4)) || gamma(end) < min(ThirdStageData(:,4))
%     disp('Third Stage Matrix Extrapolating for Trajectory Angle')
% end
% 
% if alt(end) > max(ThirdStageData(:,3)) || alt(end) < min(ThirdStageData(:,3))
%     disp('Third Stage Matrix Extrapolating for Altitude')
% end
% 
% if v(end) > max(ThirdStageData(:,5)) || v(end) < min(ThirdStageData(:,5))
%     disp('Third Stage Matrix Extrapolating for Velocity')
% end

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



[~,~,~,~,~,~, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,q11,flapdeflection,heating_rate] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,lat,lon,Alpha,eta,1, mFuel,mFuel(1),mFuel(end), 1);
[~,~,~,~,~,~, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,q12,flapdeflection2,heating_rate2] = VehicleModelCombined(gamma2, alt2, v2,auxdata,zeta2,lat2,lon2,Alpha2,eta2,throttle2, mFuel2,0,0, 0);

throttle2(M2<5.1) = 0; % remove nonsense throttle points
Isp2(M2<5.1) = 0; % remove nonsense throttle points

% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
end

% Separation_LD = lift(end)/Fd(end)

figure(2010)

subplot(5,5,[1,10])
hold on
plot(H, alt)
% plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')

dim = [.7 .52 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1562 - mFuel(end)) ' kg']},'FitBoxToText','on');  


subplot(5,5,11)
hold on
plot(time, v)

title('Velocity (m/s)')


subplot(5,5,12)
plot(time, M1)
title('Mach no')

subplot(5,5,13)
plot(time, q1)
title('Dynamic Pressure (pa)')

subplot(5,5,14)
hold on
plot(time, rad2deg(gamma))

title('Trajectory Angle (Deg)')



subplot(5,5,15)
plot(time, Fd1)
title('Drag Force')

subplot(5,5,16)
hold on
plot(time, mFuel + 8755.1 - 994)
title('Vehicle Mass (kg)')



subplot(5,5,17)
plot(time, T1)
title('Thrust (N)')

% Isp1 = T1./Fueldt1./9.81;
IspNet1 = (T1-Fd1)./Fueldt1./9.81;

subplot(5,5,18)
plot(time, Isp1)
title('Isp')

subplot(5,5,19)
plot(time, IspNet1)
title('Net Isp')

subplot(5,5,20)
plot(time, flapdeflection)
title('Flap Deflection (deg)')

subplot(5,5,21)
plot(time, rad2deg(Alpha))
title('Angle of Attack (deg)')

subplot(5,5,22)
plot(time, rad2deg(eta))
title('Bank Angle (deg)')

subplot(5,5,23)
plot(time, q11)
title('Dynamic pressure after shock')

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

line(time, M1,'Parent',ax2,'Color','k', 'LineStyle','--')

line(time, v./(10^3),'Parent',ax2,'Color','k', 'LineStyle','-.')

line(time, q1./(10^4),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

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
line(time, [rad2deg(Alpha(1:end-1)) rad2deg(Alpha(end-1))],'Parent',ax3,'Color','k', 'LineStyle','-')
line(time, [rad2deg(eta(1:end-1)) rad2deg(eta(end-1))],'Parent',ax3,'Color','k', 'LineStyle','-.')

line(time, flapdeflection,'Parent',ax3,'Color','k', 'LineStyle','--')


% line(time, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
% line(time, eq.*10,'Parent',ax3,'Color','k', 'LineStyle','-.')

line(time, IspNet1./(10^2),'Parent',ax3,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax3, 'AoA (degrees)', 'Bank Angle (degrees)','Flap Deflection (degrees)', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)


SecondStageStates = [[time time2]' [alt alt2]' [lon lon2]' [lat lat2]' [v v2]' [gamma gamma2]' [zeta zeta2]' [Alpha Alpha2]' [eta eta2]' [mFuel mFuel2]'];
dlmwrite('SecondStageStates',['time (s) ' 'altitude (m) ' 'longitude (rad) ' 'latitude (rad) ' 'velocity (m/s) ' 'trajectory angle (rad) ' 'heading angle (rad) ' 'angle of attack (rad) ' 'bank angle (rad) ' 'fuel mass (kg) '],'');
dlmwrite('SecondStageStates',SecondStageStates,'-append','delimiter',' ');
copyfile('SecondStageStates',sprintf('../ArchivedResults/%s/SecondStage_%s',Timestamp,Timestamp))



%% PLOT RETURN
addpath('..\SecondStageReturn\addaxis')

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,1)
hold on


 plot(time2,alt2/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
xlabel('time (s)');
addaxis(time2,v2/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Velocity (km/s)');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');

addaxis(time2,rad2deg(zeta2),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Heading Angle (Deg)');


legend(  'Altitude', 'Velocity', 'Heading Angle', 'location', 'best');

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,2)
hold on
plot(time2,rad2deg(Alpha2),'-','color','k', 'linewidth', 1.);
ylabel('Angle of Attack (deg)');
xlabel('time (s)')
throttle(M2<5.0)=0;
addaxis(time2,throttle2*100,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Throttle (%)');

% addaxis(time,mfuel,'-.','color','k', 'linewidth', 1.);
% addaxislabel(3,'Fuel Mass (kg)');

addaxis(time2,rad2deg(eta2),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Bank Angle (Deg)');

addaxis(time2,flapdeflection2,'-.','color','k', 'linewidth', 1.2);
addaxislabel(4,'Flap Deflection (Deg)');
legend(  'Angle of Attack', 'Throttle' , 'Bank Angle','FlapDeflection');


figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
hold on
% subplot(3,1,3)
plot(time2,M2,'-','color','k', 'linewidth', 1.);
ylabel('Mach no.')
xlabel('time (s)')

addaxis(time2,Isp2,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Specific Impulse (s)');

addaxis(time2,q2,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

% addaxis(time,L./D,':','color','k', 'linewidth', 1.);
% addaxislabel(4,'L/D');

legend(  'Mach no.', 'Isp (potential)', 'q' );


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,eta,f_t),1,mFuel(1),mFuel(end)),time(1:end),forward0);

figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,gamma);


subplot(7,1,4)
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



% Return Forward
forward0 = [alt2(1),gamma2(1),v2(1),zeta2(1),lat2(1),lon2(1), mFuel2(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time2,Alpha2,f_t),ControlInterp(time2,eta2,f_t),ThrottleInterp(time2,throttle2,f_t)),time2(1):time2(end),forward0);

% altitude  = (output.result.solution.phase(1).state(:,1)-auxdata.Re);
figure(213)
subplot(7,1,1)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time2,alt2);

% gamma  = output.result.solution.phase.state(:,5);

subplot(7,1,2)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time2,gamma2);

% latitude  = output.result.solution.phase.state(:,3);
subplot(7,1,3:5)
hold on
plot(f_y(:,6),f_y(:,5));
plot(lon2,lat2);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time2,mFuel2);

%% Check KKT and pontryagins minimum
% Check that the hamiltonian = 0 (for free end time)
% Necessary condition
input_test = output.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

lambda1 = output.result.solution.phase(1).costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test(1).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output.result.solution.phase(2).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

figure(221)
hold on
plot(time(1:end-1),H1)
plot(time2(1:end-1),H2)
ylabel('Hamiltonian')
xlabel('Time (s)')
legend('Ascent','Return')

% Check Primal Feasibility
% Check calculated derivatives with the numerical derivative of each
% porimal, scaled by that primal
figure(220)
hold on
for i = 1:length(output.result.solution.phase(1).state(1,:))
plot(time,([diff(output.result.solution.phase(1).state(:,i))./diff(output.result.solution.phase(1).time); 0] - phaseout_test(1).dynamics(:,i))./output.result.solution.phase(1).state(:,i),'--');
end
for i = 1:length(output.result.solution.phase(2).state(1,:))
    if i<= 7 % Plot different line styles when no. of colours exceeded
    plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i));
    else
    plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i),':');
    end
end
xlabel('Time (s)')
ylabel('Derivative Error')
ylim([-1,1])
legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')

%% plot engine interpolation visualiser
T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt); 
T_in1 = auxdata.interp.tempgridded(M1,rad2deg(Alpha)).*T0;
M_in1 = auxdata.interp.M1gridded(M1, rad2deg(Alpha));

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

alt3  = output.result.solution.phase(3).state(:,1);
v3    = output.result.solution.phase(3).state(:,2);
gamma3  = output.result.solution.phase(3).state(:,3);
m3    = output.result.solution.phase(3).state(:,4);
aoa3    = output.result.solution.phase(3).state(:,5);
phi3    = output.result.solution.phase(3).state(:,6);
zeta3    = output.result.solution.phase(3).state(:,7);
aoadot3       = output.result.solution.phase(3).control(:,1);

forward0 = [alt3(1),v3(1),gamma3(1),m3(1),phi3(1),zeta3(1)];

time3 = output.result.solution.phase(3).time;

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel3_forward(f_t, f_y,auxdata,ControlInterp(time3,aoa3,f_t),ControlInterp(time3,aoadot3,f_t)),time3(1:end),forward0);

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T, L, D, q] = ThirdStageDyn(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,phi3,zeta3);


xi(1) = lon(end);
for i = 2:length(time3)
    xi(i) = xi(i-1) + xidot(i-1)*(time3(i)-time3(i-1));
end

[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,phiexo,incexo,Texo,CLexo,Lexo,incdiffexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), phi3(end),xi(end), zeta3(end), m3(end), auxdata);




figure(312)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time3,alt3);

figure(313)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time3,v3);


figure(301)
    addpath('addaxis')
    hold on

    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [alt3; altexo.']/1000, 'LineStyle', '-','Color','k', 'lineWidth', 2.2)
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[q;qexo.';qexo(end)]/1000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(aoa3);0*ones(length(timeexo),1)], 'LineStyle', '--','Color','k', 'lineWidth', 0.7)
    ylabel('Altitude (km), Dynamic Pressure (kPa), Angle of Attack (deg)');
    

    addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[v3;v3exo.'], [0 7000], 'LineStyle', '--','Color','k', 'lineWidth', 1.8)
    addaxisplot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[ m3;mexo.';mexo(end)],2, 'LineStyle', ':','Color','k', 'lineWidth', 1.3)
    addaxislabel(2,'Velocity (m/s), Mass (kg)');


    addaxis([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(Vec_angle);0*ones(length(timeexo),1)], 'LineStyle', ':','Color','k', 'lineWidth', 2.1)
    addaxisplot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [rad2deg(gamma3);rad2deg(gammaexo).'],3, 'LineStyle', '-','Color','k', 'lineWidth', .6)
    addaxislabel(3,'Thrust Vector Angle (deg), Trajectory Angle (deg)');

    legend(  'Altitude','Dynamic Pressure','Angle of Attack', 'Velocity',  'Mass', 'Thrust Vector Angle', 'Trajectory Angle' );
    xlabel('Time (s)');
    xlim([0 timeexo(end)+time3(end)-time3(1)])
    box off
    % Write data to file
    dlmwrite('ThirdStageData',[time3, alt3, v3, m3,q ,gamma3,D ,zeta3], ' ')

    
%% SAVE FIGS

saveas(figure(301),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ThirdStage.fig']);
print(figure(301),'ThirdStage','-dpng');
movefile('ThirdStage.png',sprintf('../ArchivedResults/%s',Timestamp));
saveas(figure(202),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'SecondStage.fig']);
saveas(figure(1),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Return1.fig']);    
saveas(figure(2),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Return2.fig']);
saveas(figure(3),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Return3.fig']);
saveas(figure(221),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Hamiltonian.fig']);
saveas(figure(220),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Validation.fig']);
saveas(figure(212),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Forward1.fig']);
saveas(figure(213),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'Forward2.fig']);
% saveas(figure(2100),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'eq.fig']);
% saveas(figure(2110),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ISP.fig']);
saveas(figure(2301),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'GroundTrack.fig']);




%% Run First Stage =========================================================
% const_firststage = 1;
% addpath('../../FirstStage')
% % addpath('../../DIDO_7.3.7')
% % run startup.m
% [FirstStageStates] = FirstStageProblem(alt(1),gamma(1),lat(1),zeta(1),const_firststage);
% % cd('../SecondStage/Combined - 2nd Stage Ascent and Return')
% dlmwrite('FirstStageStates', FirstStageStates);
% copyfile('FirstStageStates',sprintf('../ArchivedResults/%s/firststage_%s',Timestamp,Timestamp))
% % Latitude Plot
% figure(250)
% plot(FirstStageStates(:,9))
% plot(phi)
% plot(ThirdStagePhi)
% title('Latitude')
% 
% saveas(figure(101),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'FirstStage.fig']);
%% Create Easy Latex Inputs



sprintf('\newcommand{\FirstStageAltMode', mode ,'}{ ', FirstStageStates(1,end) , '}');



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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket-Scramjet-Rocket Launch Optimiser
% By Sholto Forbes-Spyratos
% Utilises the GPOPS-2 proprietary optimisation software
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

mode = 2
% auxdata.mode = mode;
returnMode = 1

%% Misc Modifiers

auxdata.delta = deg2rad(0) % thrust vector angle test
auxdata.dragmod = 1. %drag increase test

%% Launch Point
lat0 = deg2rad(-12.164); % Equatorial Launch Australia Spaceport near Nhulunbuy
lon0 = deg2rad(136.755);



%% FIRST STAGE


% auxdata.Throttle = 0.75; % throttle the Merlin engine down by a modeant value, to enable easier pitchover
 auxdata.Throttle = .7
% Aerodynamics File Path
Aero = dlmread('FirstStageAeroCoeffs.txt');


%% Vehicle 


% mRocket =21816; % total mass of scaled Falcon at 9.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx

mRocket =19569; % total mass of scaled Falcon at 8.5m, note, this will not be the final total mass. Calculated using the method outlined in SIZING.docx

mEngine = 470; % Mass of Merlin 1C
mFuel = 0.939*(mRocket-mEngine); % structural mass fraction calculated without engine
mSpartan = 9819.11;

% Thrust and Isp are modified with altitude through the formula:
% SL + (101325-P_atm)*Mod

auxdata.Vehicle.T.SL = 555900; % Thrust from Falcon 1 users guide. 
auxdata.Vehicle.T.Mod = 0.5518; % exit area calculated in SIZING.docx

auxdata.Vehicle.Isp.SL = 275; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
auxdata.Vehicle.Isp.Mod = 2.9410e-04;

auxdata.Vehicle.Area = 62.77; 


%% Import Atmosphere

auxdata.Atmosphere = dlmread('atmosphere.txt');

%% Calculate Aerodynamic Splines

interp.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
interp.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = interp.Lift(grid.M,grid.AoA);
grid.Drag = interp.Drag(grid.M,grid.AoA);

auxdata.interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
auxdata.interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');

mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)


%% Assign Pitchover Conditions

%Define initial conditions at pitchover, these are assumed
h0 = 90;  
v0 = 30;    

gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

mF = mEmpty+mSpartan;  %Assume that we use all of the fuel


alpha0 = 0; %Set initial angle of attack to 0

hMin1 = 1;   %Cannot go through the earth
hMax1 = 30000;  

vMin1 = 0; 
vMax1 = 3000;  

mMin1 = mEmpty;
mMax1 = mTotal;

phiMin1 = -0.5;
phiMax1 = 0.5;

zetaMin1 = -2*pi;
zetaMax1 = 2*pi;

alphaMin1 = -deg2rad(5);
alphaMax1 = deg2rad(2);

dalphadt2Min1 = -0.1;
dalphadt2Max1 = 0.1;

lonMin = 2;         lonMax = 3;

gammaMin1 = deg2rad(-.1);
gammaMax1 = gamma0;
% This sets the control limits, this is second derivative of AoA
uMin1 = [-.0005]; % Can do either AoA or thrust
uMax1 = [.0005];

%-------------------------------------------
% Set up the problem bounds in SCALED units
%-------------------------------------------

tfMax1 	    = 300;     % large upper bound; do not choose Inf
		 
%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%

bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;

bounds.phase(1).finaltime.lower = 0;
bounds.phase(1).finaltime.upper = tfMax1;

bounds.phase(1).initialstate.lower = [h0, v0,  mF-1, gamma0, alpha0,  zetaMin1, dalphadt2Min1, lat0, lon0 ];
bounds.phase(1).initialstate.upper = [h0, v0, mMax1, gamma0, alpha0, zetaMax1, dalphadt2Max1, lat0, lon0];

bounds.phase(1).state.lower = [hMin1, vMin1, mF-1, gammaMin1, alphaMin1, zetaMin1, dalphadt2Min1, phiMin1, lonMin ];
bounds.phase(1).state.upper = [ hMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax1, dalphadt2Max1, phiMax1, lonMax];

bounds.phase(1).finalstate.lower = [hMin1, vMin1, mF-1, gammaMin1, alphaMin1, zetaMin1, dalphadt2Min1, phiMin1, lonMin ];
bounds.phase(1).finalstate.upper = [ hMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax1, dalphadt2Max1, phiMax1, lonMax];

bounds.phase(1).control.lower = uMin1;
bounds.phase(1).control.upper = uMax1;

bounds.phase(1).path.lower = 0;
bounds.phase(1).path.upper = 50000;

% Tie stages together
bounds.eventgroup(1).lower = [zeros(1,8)];
bounds.eventgroup(1).upper = [zeros(1,8)]; 
% bounds.eventgroup(1).lower = [ones(1,7)*-100000 0];
% bounds.eventgroup(1).upper = [ones(1,7)*100000 0]; 

guess.phase(1).time =  [0; tfMax1];

guess.phase(1).state(:,1) = [h0; h0];
guess.phase(1).state(:,2) = [v0; 1500];
guess.phase(1).state(:,3) = [mMax1; mF];
guess.phase(1).state(:,4) = [gamma0; 0];
guess.phase(1).state(:,5) = [alpha0; 0];
guess.phase(1).state(:,6) = [0; 0];
guess.phase(1).state(:,7) = [0; 0];
guess.phase(1).state(:,8) = [lat0; lat0];
guess.phase(1).state(:,9) = [lon0; lon0];

guess.phase(1).control = [0; 0];


%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\..\thirdStage')
addpath('..\EngineData')
addpath('..\')

% Timestamp = datestr(now,30)
% mkdir('../ArchivedResults', strcat(Timestamp, 'mode', num2str(mode)))
% copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',strcat(Timestamp,'mode',num2str(mode))))

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

aero2.aero_EngineOff = importdata('SPARTANaero15.1557');
aero2.flapaero = importdata('SPARTANaeroFlaps15.1557');
aero2.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1557');
aero2.aero_Engine = importdata('SPARTANEngine15.1557');
aero2.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero2.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.cylTankEnd,auxdata.interp.Cd_spline_EngineOff.cylTankEnd,auxdata.interp.Cl_spline_EngineOn.cylTankEnd,auxdata.interp.Cd_spline_EngineOn.cylTankEnd,auxdata.interp.flap_spline_EngineOff.cylTankEnd,auxdata.interp.flap_spline_EngineOn.cylTankEnd] = AeroInt(aero2,auxdata,T_L,CG_z);



% End of acceleration, with third stage. 
% NOTE: it is assumed that the CG does not change
% due to fuel after the end of acceleration phase, so there will still be some fuel left when this is used. 

CG_z = (-0.2134*4.9571e+03+ 3300*0.547)/(4.9571e+03+3300);

aero3.aero_EngineOff = importdata('SPARTANaero15.727');
aero3.flapaero = importdata('SPARTANaeroFlaps15.727');
aero3.aero_EngineOn = importdata('SPARTANaeroEngineOn15.727');
aero3.aero_Engine = importdata('SPARTANEngine15.727');
aero3.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero3.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero3,auxdata,T_L,CG_z);

% Flyback, without third stage. Fuel variation not used for flyback.

CG_z = -0.2134; % calculated fom CREO

aero4.aero_EngineOff = importdata('SPARTANaero15.1255');
aero4.flapaero = importdata('SPARTANaeroFlaps15.1255');
aero4.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1255');
aero4.aero_Engine = importdata('SPARTANEngine15.1255');
aero4.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero4.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStage,auxdata.interp.Cd_spline_EngineOff.noThirdStage,auxdata.interp.Cl_spline_EngineOn.noThirdStage,auxdata.interp.Cd_spline_EngineOn.noThirdStage,auxdata.interp.flap_spline_EngineOff.noThirdStage,auxdata.interp.flap_spline_EngineOn.noThirdStage] = AeroInt(aero4,auxdata,T_L,CG_z);


%% Import Bounds %%========================================================

latMin2 = -0.5;  latMax2 = 0.5;

% lat0 = -0.264;
% lon0 = deg2rad(145);

aoaMin21 = 0;  aoaMax21 = 10*pi/180;
if returnMode == 0
bankMin21 = -1*pi/180; bankMax21 =   1*pi/180;    
else
bankMin21 = -1*pi/180; bankMax21 =   90*pi/180;
end
% Primal Bounds
bounds.phase(2).state.lower = [Stage2.Bounds.Alt(1), lonMin, latMin2, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Bounds.mFuel(1)];
bounds.phase(2).state.upper = [Stage2.Bounds.Alt(2), lonMax, latMax2, Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Bounds.mFuel(2)];

% Initial States
bounds.phase(2).initialstate.lower = [Stage2.Bounds.Alt(1),lonMin, latMin2, Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
bounds.phase(2).initialstate.upper = [Stage2.Bounds.Alt(2),lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(15), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Initial.mFuel];
% bounds.phase(2).initialstate.lower = [Stage2.Bounds.Alt(1),lon0, lat0, 1500, Stage2.Bounds.gamma(1), Stage2.Bounds.zeta(1), aoaMin21, bankMin21, Stage2.Initial.mFuel] ;
% bounds.phase(2).initialstate.upper = [Stage2.Bounds.Alt(2),lon0, lat0, 1500, deg2rad(15), Stage2.Bounds.zeta(2), aoaMax21, bankMax21, Stage2.Initial.mFuel];

% End States
% End bounds are set slightly differently, to encourage an optimal solution
bounds.phase(2).finalstate.lower = [20000, lonMin, latMin2, 2300, 0, Stage2.Bounds.zeta(1), aoaMin21, 0, Stage2.End.mFuel];
bounds.phase(2).finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), deg2rad(20), Stage2.Bounds.zeta(2), aoaMax21, 0, Stage2.Initial.mFuel];

% bounds.phase(2).finalstate.lower = [34000, lonMin, latMin2, 2300, 0, Stage2.Bounds.zeta(1), aoaMin21, 0, Stage2.End.mFuel];
% bounds.phase(2).finalstate.upper = [45000, lonMax, latMax2, Stage2.Bounds.v(2), 0, Stage2.Bounds.zeta(2), aoaMax21, 0, Stage2.Initial.mFuel];
%  disp('gamma set to 0')
 
% Control Bounds
bounds.phase(2).control.lower = [deg2rad(-.2), deg2rad(-1)];
bounds.phase(2).control.upper = [deg2rad(.2), deg2rad(1)];

% Time Bounds
bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = Stage2.Bounds.time(2);
bounds.phase(2).finaltime.lower = Stage2.Bounds.time(1);
bounds.phase(2).finaltime.upper = Stage2.Bounds.time(2);


%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.
if mode == 1 || mode == 14 || mode == 15 || mode == 2
    bounds.phase(2).path.lower = [0];
    bounds.phase(2).path.upper = [50000];

elseif mode ==3 || mode == 32
        bounds.phase(2).path.lower = [49970];
    bounds.phase(2).path.upper = [50030];
end

% bounds.phase(2).integral.lower = -10000000000;
% bounds.phase(2).integral.upper = 10000000000;
%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
% guess.phase(2).state(:,1)   = [24000;35000];
guess.phase(2).state(:,2)   = [2.53;2.5368];
guess.phase(2).state(:,3)   = [-0.269;-0.10];
guess.phase(2).state(:,4)   = Stage2.Guess.v.';
guess.phase(2).state(:,5)   = Stage2.Guess.gamma.';
guess.phase(2).state(:,6)   = Stage2.Guess.zeta.';
guess.phase(2).state(:,7)   = [2*pi/180; 5*pi/180];
guess.phase(2).state(:,8)   = [deg2rad(10);deg2rad(10)];
guess.phase(2).state(:,9) 	= [Stage2.Initial.mFuel; 100];

guess.phase(2).control      = [[0;0],[0;0]];
guess.phase(2).time          = [0;650];

% guess.phase(2).integral = 0

% Tie stages together
bounds.eventgroup(2).lower = [zeros(1,9)];
bounds.eventgroup(2).upper = [zeros(1,9)]; 

%% Flyback
tfMin = 0;            tfMax = 5000;
altMin = 10;  altMax = 70000;
speedMin = 10;        speedMax = 5000;
fpaMin = -80*pi/180;  fpaMax =  80*pi/180;
aziMin = 60*pi/180; aziMax =  500*pi/180;
mFuelMin = 0; mFuelMax = 500;
bankMin21 = -10*pi/180; bankMax21 =   90*pi/180;

throttleMin = 0; throttleMax = 1;

bounds.phase(3).initialtime.lower = 0;
bounds.phase(3).initialtime.upper = 3000;
bounds.phase(3).finaltime.lower = 400;
bounds.phase(3).finaltime.upper = 4000;

% Initial Bounds

 bounds.phase(3).initialstate.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, 0, mFuelMin, throttleMin];
bounds.phase(3).initialstate.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, 0, mFuelMax, throttleMax];   

% State Bounds
if returnMode == 0
bounds.phase(3).state.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, 1, throttleMax];
else
bounds.phase(3).state.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
end
% End State Bounds
if returnMode == 0
bounds.phase(3).finalstate.lower = [altMin, lonMin, latMin2, speedMin, fpaMin, aziMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).finalstate.upper = [altMax, lonMax, latMax2, speedMax, fpaMax, aziMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
else
bounds.phase(3).finalstate.lower = [altMin, lon0, lat0, speedMin, deg2rad(-20), aziMin, aoaMin21, bankMin21, Stage2.End.mFuel, throttleMin];
bounds.phase(3).finalstate.upper = [altMax, lon0, lat0, speedMax, deg2rad(30), aziMax, aoaMax21, bankMax21, Stage2.End.mFuel, throttleMax];
end
% Control Bounds
bounds.phase(3).control.lower = [deg2rad(-.2), deg2rad(-5), -1];
bounds.phase(3).control.upper = [deg2rad(.2), deg2rad(5), 1];

% Path Bounds
bounds.phase(3).path.lower = 0;
bounds.phase(3).path.upper = 50000;

bounds.eventgroup(3).lower = [0]; % Constrain final latitude and longitude, with variable first stage approximation
bounds.eventgroup(3).upper = [0]; 

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
mFuelGuess          = [100; mFuelMin];
guess.phase(3).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[0.;0.]];
guess.phase(3).control = [[0;0],[0;0],[0;0]];
guess.phase(3).time    = tGuess;


%% Third Stage

% tfMin3 = 7;            tfMax3 = 200;
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
bounds.phase(4).initialtime.lower = 0;
bounds.phase(4).initialtime.upper = 10000;
bounds.phase(4).finaltime.lower = 1;
bounds.phase(4).finaltime.upper = 10000;


bounds.phase(4).initialstate.lower = [altMin3,vMin3, 0,  auxdata.Stage3.mTot, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).initialstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(4).state.lower = [altMin3,vMin3, gammaMin3, 0, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).state.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

bounds.phase(4).finalstate.lower = [altMin3, vMin3, 0, 0, 0, phiMin3, zetaMin3];
bounds.phase(4).finalstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, 0, phiMax3, zetaMax3];

bounds.phase(4).control.lower = [aoadotMin3];
bounds.phase(4).control.upper = [aoadotMax3];

bounds.phase(4).path.lower = [-deg2rad(8), -inf];
bounds.phase(4).path.upper = [deg2rad(8), 0];


bounds.eventgroup(4).lower = [zeros(1,7) 90000 0];
bounds.eventgroup(4).upper = [zeros(1,6) 1000 566000 0];

tGuess              = [0; 150];
altGuess            = [35000; 60000];
vGuess          = [2700; 6000];
gammaGuess            = [0; deg2rad(10)];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(20); deg2rad(20)];
phiGuess = [-0.11;-0.11];
zetaGuess = [deg2rad(97);deg2rad(97)];
guess.phase(4).state   = [altGuess, vGuess, gammaGuess, mGuess, aoaGuess, phiGuess, zetaGuess];
guess.phase(4).control = [0;0];
guess.phase(4).time    = tGuess;

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
% mesh.method       = 'hp-LiuRao-Legendre'; % Default method does not perform h adaptions sometimes, and can be unstable. 
% mesh.method       = 'hp-LiuRao'; % Default method does not perform h adaptions sometimes, and can be unstable. 
 mesh.method       = 'hp-DarbyRao';
mesh.maxiterations = 1;
mesh.colpointsmin = 8;
mesh.colpointsmax = 50;
mesh.tolerance    = 1e-5;

% mesh.phase(2).fraction = 0.1*ones(1,10)
% 
% mesh.phase(3).fraction = 0.05*ones(1,20)
% 
% mesh.phase(4).fraction = 0.1*ones(1,10)
% 
% mesh.phase(2).colpoints = 4*ones(1,10)
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

setup.nlp.ipoptoptions.maxiterations = 100;

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



if mode == 1
    setup_variations{1} = setup;
    
elseif mode == 2
    q_vars = [40000 50000]
    for i = 1:length(q_vars)  

    setup_variations{i} = setup;
    setup_variations{i}.bounds.phase(2).path.upper = q_vars(i);
    setup_variations{i}.bounds.phase(3).path.upper = q_vars(i);
    end
    
end



for j = 1:length(setup_variations)

    
for i = 1:num_it
setup_par(i) = setup_variations{j};
% setup_par(i).nlp.ipoptoptions.maxiterations = 500 + 10*i;
setup_par(i).guess.phase(2).state(:,1)   = [24000;24000 + 1000*i]; % vary altitude guess
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

PayloadMass(i) = -output_temp.result.objective;

output_store{i} = output_temp;

end

% [min_error,index] = min(error); % Calculate the result which minimises the chosen error function

[max_pl,index] = max(PayloadMass);% Calculate the result which maximises payload mass the chosen error function
output{j} = output_store{index};

end

%%

% EndTime = datestr(now,30) % Display the ending time

SinglePlotter(output{1},auxdata,mode,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,mRocket,mSpartan,mFuel,h0,v0,bounds);



% ----------- Reusable Launch Vehicle Entry Example ------------%
% This example is taken verbatim from the following reference:  %
% Betts, J. T., Practical Methods for Optimal Control Using     %
% Nonlinear Programming, SIAM Press, Philadelphia, 2009.        %
% --------------------------------------------------------------%
%close all
clear all
clc
addpath('..\')
addpath('..\EngineData')



auxdata.const = 1 % 1 normal, 2 10% increased Cd, 3 10% decreased Cd, 4 10% increased Isp, 5: 10% decreased Isp, 6:10% increased Cl, 7:10% decreased Cl, 8:10% increased range

%% Inputs ============================================

%% Aerodynamic Data - Engine off


% flapaero = importdata('SPARTAN_Flaps.txt');
% 
% interp.flap_momentCl_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,3), 'linear', 'nearest');
% interp.flap_momentCd_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,4), 'linear', 'nearest');
% interp.flap_momentdef_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,2), 'linear', 'nearest');

aero_EngineOff = importdata('SPARTANaero-14-5mCG.txt');

interp.Cl_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,3));
interp.Cd_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,4));
interp.Cm_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,5));

[MList_EngineOff,AOAList_EngineOff] = ndgrid(unique(aero_EngineOff(:,1)),unique(aero_EngineOff(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid_EngineOff = [];
Cd_Grid_EngineOff = [];
Cm_Grid_EngineOff = [];
flap_Grid = [];

for i = 1:numel(MList_EngineOff)
    M_temp = MList_EngineOff(i);
    AoA_temp = AOAList_EngineOff(i);
    
    Cl_temp = interp.Cl_scattered_EngineOff(M_temp,AoA_temp);
    Cd_temp = interp.Cd_scattered_EngineOff(M_temp,AoA_temp);
    Cm_temp = interp.Cm_scattered_EngineOff(M_temp,AoA_temp);
    
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOff(M_temp,0);
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOff(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOff(M_temp,0);
    
%     Cl_AoA0_withflaps_temp = interp.flap_momentCl_scattered(M_temp,-(Cm_temp-Cm_temp_AoA0));
%     Cd_AoA0_withflaps_temp = interp.flap_momentCd_scattered(M_temp,-(Cm_temp-Cm_temp_AoA0)) ;
%     
%     flap_Cl_temp = Cl_AoA0_withflaps_temp - Cl_temp_AoA0;
%     flap_Cd_temp = Cd_AoA0_withflaps_temp - Cd_temp_AoA0;
    
    I = cell(1, ndims(MList_EngineOff)); 
    [I{:}] = ind2sub(size(MList_EngineOff),i);
    
%     Cl_Grid(I{(1)},I{(2)}) = Cl_temp+flap_Cl_temp;
%     Cd_Grid(I{(1)},I{(2)}) = Cd_temp+flap_Cd_temp;
%     Cm_Grid(I{(1)},I{(2)}) = Cm_temp;
% 
%     flap_Grid(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp-Cm_temp_AoA0)) ;
%     
%     Cl_Grid_test(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_test(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_test(I{(1)},I{(2)}) = Cm_temp;
    
    Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp;
    Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp;
    Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp;
end
auxdata.interp.Cl_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cl_Grid_EngineOff,'spline','linear');
auxdata.interp.Cd_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cd_Grid_EngineOff,'spline','linear');
auxdata.interp.Cm_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cm_Grid_EngineOff,'spline','nearest');

%% Aerodynamic Data - Engine on 

% flapaero_154 = importdata('SPARTAN_Flaps-15-4mCG.txt');
% 
% interp.flap_momentCl_scattered_154 = scatteredInterpolant(flapaero_154(:,1),flapaero_154(:,5),flapaero_154(:,3), 'linear', 'nearest'); % interpolation by moments
% interp.flap_momentCd_scattered_154 = scatteredInterpolant(flapaero_154(:,1),flapaero_154(:,5),flapaero_154(:,4), 'linear', 'nearest');
% interp.flap_momentdef_scattered_154 = scatteredInterpolant(flapaero_154(:,1),flapaero_154(:,5),flapaero_154(:,2), 'linear', 'nearest');

aero_EngineOn = importdata('SPARTANaero-EngineOn-15-4mCG.txt');


interp.Cl_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3));
interp.Cd_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,4));
interp.Cm_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,5));

[MList_EngineOn,AOAList_EngineOn] = ndgrid(unique(aero_EngineOn(:,1)),unique(aero_EngineOn(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid_EngineOn = [];
Cd_Grid_EngineOn = [];
Cm_Grid_EngineOn = [];
flap_Grid = [];

for i = 1:numel(MList_EngineOn)
    M_temp = MList_EngineOn(i);
    AoA_temp = AOAList_EngineOn(i);
    
    Cl_temp = interp.Cl_scattered_EngineOn(M_temp,AoA_temp);
    Cd_temp = interp.Cd_scattered_EngineOn(M_temp,AoA_temp);
    Cm_temp = interp.Cm_scattered_EngineOn(M_temp,AoA_temp);
    
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOn(M_temp,0);
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOn(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOn(M_temp,0);
    
%     Cl_AoA0_withflaps_temp = interp.flap_momentCl_scattered_154(M_temp,-(Cm_temp-Cm_temp_AoA0));
%     Cd_AoA0_withflaps_temp = interp.flap_momentCd_scattered_154(M_temp,-(Cm_temp-Cm_temp_AoA0)) ;
%     
%     flap_Cl_temp = Cl_AoA0_withflaps_temp - Cl_temp_AoA0;
%     flap_Cd_temp = Cd_AoA0_withflaps_temp - Cd_temp_AoA0;
    
    I = cell(1, ndims(MList_EngineOn)); 
    [I{:}] = ind2sub(size(MList_EngineOn),i);
    
%     Cl_Grid(I{(1)},I{(2)}) = Cl_temp+flap_Cl_temp;
%     Cd_Grid(I{(1)},I{(2)}) = Cd_temp+flap_Cd_temp;
%     Cm_Grid(I{(1)},I{(2)}) = Cm_temp;
% 
%     flap_Grid(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp-Cm_temp_AoA0)) ;
%     
%     Cl_Grid_test(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_test(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_test(I{(1)},I{(2)}) = Cm_temp;
    
    Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp;
    Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp;
    Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp;
end
auxdata.interp.Cl_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cl_Grid_EngineOn,'spline','linear');
auxdata.interp.Cd_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cd_Grid_EngineOn,'spline','linear');
auxdata.interp.Cm_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cm_Grid_EngineOn,'spline','nearest');


% Produce Atmosphere Data
Atmosphere = dlmread('atmosphere.txt');
auxdata.Atmosphere = Atmosphere; 
auxdata.interp.c_spline = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data

auxdata.interp.rho_spline = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,4)); % Calculate density using atmospheric data

auxdata.interp.T0_spline = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,2)); 

auxdata.interp.P0_spline = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,3)); 




%% Engine Data
% Import engine data
auxdata.engine_data = dlmread('ENGINEDATA.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.engine_data;


%%
 p=polyfitn([engine_data(:,1),engine_data(:,2)],engine_data(:,3),3);

M_englist = unique(sort(engine_data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(engine_data(:,1)));

T_englist = unique(sort(engine_data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = unique(sort(engine_data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);
%%
% This sets the equivalence ratio interpolation region. VERY IMPORTANT

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

auxdata.equivalence = scatteredInterpolant(eq_data(:,1),eq_data(:,2),eq_data(:,4), 'linear');
grid.eq_eng = auxdata.equivalence(grid.Mgrid_eng,grid.T_eng);
auxdata.interp.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');


%% Load the interpolated Isp data

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

%% Shock Data
% Import conical shock data and create interpolation splines 
shockdata = dlmread('ShockMat');
[MList2,AOAList2] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
auxdata.interp.M1gridded = griddedInterpolant(MList2,AOAList2,M1_Grid,'spline','linear');
auxdata.interp.presgridded = griddedInterpolant(MList2,AOAList2,pres_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList2,AOAList2,temp_Grid,'spline','linear');


%% Import Vehicle and trajectory Config Data %%============================
addpath('../')
run VehicleConfig.m
run TrajectoryConfig50kPa.m

auxdata.Stage3 = Stage3;
auxdata.Stage2 = Stage2;

%=============================================== 

%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)
% auxdata.mass = mstruct;               % Vehicle Mass (kg)
% auxdata.A = 62.77; %m^2
%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0     = 0;
alt0   = 34500;   
rad0   = alt0+auxdata.Re;
altf   = 500;   
radf   = altf+auxdata.Re;


i=1
% for lon0 = deg2rad(144):deg2rad(.5):deg2rad(147)
%     for lat0 = deg2rad(-10):deg2rad(.5):deg2rad(-7.5)
%         for j = 1:4
        
%   lon0   =   2.552;    
lon0   = deg2rad(145);%ref
lonf = deg2rad(145);
% lat0 = -0.1571;
lat0   = -0.135;%ref
latf   = -0.269;

if auxdata.const ==8
    latf = latf + (latf-lat0)*0.1;
end
if auxdata.const ==9
    latf = latf - (latf-lat0)*0.1;
end

speed0 = +2878;
% speed0 = 1000;

speedf = 100;
fpa0   = 2.9*pi/180; 
fpaf   = 0*pi/180;
azi0   = +102*pi/180; 
% azi0   = +180*pi/180; 
azif   = 270*pi/180;

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
tfMin = 0;            tfMax = 5000;
radMin = auxdata.Re;  radMax = auxdata.Re+70000;
lonMin = -pi;         lonMax = -lonMin;
latMin = -70*pi/180;  latMax = -latMin;
speedMin = 10;        speedMax = 5000;
fpaMin = -80*pi/180;  fpaMax =  80*pi/180;
aziMin = 60*pi/180; aziMax =  360*pi/180;
mFuelMin = 0; mFuelMax = 300;

aoaMin = 0;  aoaMax = 10*pi/180;
bankMin = -1*pi/180; bankMax =   90*pi/180;
throttleMin = 0; throttleMax = 1;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
% bounds.phase.initialtime.lower = t0;
% bounds.phase.initialtime.upper = t0;
% bounds.phase.finaltime.lower = tfMin;
% bounds.phase.finaltime.upper = tfMax;
% bounds.phase.initialstate.lower = [rad0, lon0, lat0, speed0, fpa0, azi0];
% bounds.phase.initialstate.upper = [rad0, lon0, lat0, speed0, fpa0, azi0];
% bounds.phase.state.lower = [radMin, lonMin, latMin, speedMin, fpaMin, aziMin];
% bounds.phase.state.upper = [radMax, lonMax, latMax, speedMax, fpaMax, aziMax];
% bounds.phase.finalstate.lower = [radMin, lonMin, latMin, speedMin, deg2rad(-5), azif];
% bounds.phase.finalstate.upper = [radf, lonMax, latMax, speedMax, deg2rad(5), azif];
% bounds.phase.control.lower = [aoaMin, bankMin];
% bounds.phase.control.upper = [aoaMax, bankMax];

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [rad0, lon0, lat0, speed0, fpa0, azi0, aoaMin, bankMin, mFuelMin, throttleMin];
bounds.phase.initialstate.upper = [rad0, lon0, lat0, speed0, fpa0, azi0, aoaMax, bankMax, mFuelMax, throttleMax];

bounds.phase.state.lower = [radMin, lonMin, latMin, speedMin, fpaMin, aziMin, aoaMin, bankMin, mFuelMin, throttleMin];
bounds.phase.state.upper = [radMax, lonMax, latMax, speedMax, fpaMax, aziMax, aoaMax, bankMax, mFuelMax, throttleMax];

bounds.phase.finalstate.lower = [radMin, lonf-0.002, latf-0.002, speedMin, deg2rad(-10), aziMin, aoaMin, bankMin, mFuelMin, throttleMin];
bounds.phase.finalstate.upper = [200+auxdata.Re, lonf+0.002, latf+0.002, speedMax, deg2rad(30), aziMax, aoaMax, bankMax, mFuelMin, throttleMax];
% bounds.phase.finalstate.lower = [radMin, lonf-0.002, latf-0.002, speedMin, deg2rad(-60), aziMin, aoaMin, bankMin, mFuelMin, throttleMin];
% bounds.phase.finalstate.upper = [200+auxdata.Re, lonf+0.002, latf+0.002, speedMax, deg2rad(60), aziMax, aoaMax, bankMax, mFuelMin, throttleMax];

% bounds.phase.finalstate.lower = [radMin,lonf-0.002, latMin, speedMin, deg2rad(-10), deg2rad(260), aoaMin, bankMin, mFuelMin, throttleMin];
% bounds.phase.finalstate.upper = [radMax, lonMax, latMax, speedMax, deg2rad(30), deg2rad(280), aoaMax, bankMax, mFuelMin, throttleMax];
% bounds.phase.finalstate.lower = [radMin, lonf-0.002, latMin, speedMin, deg2rad(-10), aziMin, aoaMin, bankMin, mFuelMin, throttleMin];
% bounds.phase.finalstate.upper = [200+auxdata.Re, lonf+0.002, latMax, speedMax, deg2rad(60),  aziMax, aoaMax, bankMax, mFuelMin, throttleMax];

bounds.phase.control.lower = [deg2rad(-.5), deg2rad(-5), -1];
bounds.phase.control.upper = [deg2rad(.5), deg2rad(5), 1];

bounds.phase.path.lower = 0;
bounds.phase.path.upper = 50000;

% bounds.phase.path.lower = [0, 0];
% bounds.phase.path.upper = [50000, 0];

% bounds.phase.path.lower = [0, -15];
% bounds.phase.path.upper = [50000, 15];
%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
% tGuess              = [0; 500];
% radGuess            = [rad0; radf];
% lonGuess            = [lon0; lon0+10*pi/180];
% latGuess            = [lat0; lat0+10*pi/180];
% speedGuess          = [speed0; speedf];
% fpaGuess            = [fpa0; fpaf];
% aziGuess            = [azi0; azif];
% aoaGuess            = [8*pi/180; 8*pi/180];
% bankGuess           = [60*pi/180; 60*pi/180];
% guess.phase.state   = [radGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess];
% guess.phase.control = [aoaGuess, bankGuess];
% % guess.phase.control = [aoaGuess];
% guess.phase.time    = tGuess;

tGuess              = [0; 1000];
radGuess            = [rad0; radf];
lonGuess            = [lon0; lon0-.1*pi/180];

if auxdata.const == 4
latGuess            = [lat0; lat0];
else
latGuess            = [lat0; lat0+0.1*pi/180];
end

speedGuess          = [speed0; speedf];
fpaGuess            = [fpa0; fpaf];
aziGuess            = [azi0; azif];

if auxdata.const ==1
aoaGuess            = [6*pi/180; 6*pi/180];
bankGuess           = [75*pi/180; 0*pi/180];
elseif auxdata.const == 2

    aoaGuess            = [3*pi/180; 3*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];
elseif auxdata.const == 3
     aoaGuess            = [6*pi/180; 6*pi/180];
bankGuess           = [89*pi/180; 0*pi/180]; 
elseif auxdata.const == 4
     aoaGuess            = [7*pi/180; 3*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];  

elseif auxdata.const == 5
     aoaGuess            = [8*pi/180; 2*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];  
elseif auxdata.const == 6
     aoaGuess            = [8*pi/180; 2*pi/180];
bankGuess           = [89*pi/180; 80*pi/180];
elseif auxdata.const == 7
     aoaGuess            = [5*pi/180; 2*pi/180];
bankGuess           = [89*pi/180; 80*pi/180];
elseif auxdata.const == 8
     aoaGuess            = [3*pi/180; 3*pi/180];
bankGuess           = [89*pi/180; 0*pi/180];
elseif auxdata.const == 9
     aoaGuess            = [5*pi/180; 5*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];
end
% mFuelGuess          = [mFuelMax; mFuelMin];
mFuelGuess          = [100; mFuelMin];
guess.phase.state   = [radGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[0;0]];
guess.phase.control = [[0;0],[0;0],[0;0]];
% guess.phase.control = [aoaGuess];
guess.phase.time    = tGuess;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
% if auxdata.const == 9 
% mesh.maxiterations = 8;
% else
    mesh.maxiterations = 6;
%     end
mesh.colpointsmin = 2;
mesh.colpointsmax = 15;
mesh.tolerance    = 1e-5;


%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @SecondStageReturnContinuous;
setup.functions.endpoint             = @SecondStageReturnEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
% setup.derivatives.supplier = 'adigator';
% setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
if auxdata.const == 3 || auxdata.const == 6 || auxdata.const == 7 || auxdata.const == 9 
   setup.nlp.ipoptoptions.maxiterations = 400; 
% elseif auxdata.const == 9 
%     setup.nlp.ipoptoptions.maxiterations = 350; 
else
setup.nlp.ipoptoptions.maxiterations = 300;
end
%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);


solution = output.result.solution;
rad  = solution.phase.state(:,1);
lon  = solution.phase.state(:,2);
lat  = solution.phase.state(:,3);
v    = solution.phase.state(:,4);
fpa  = solution.phase.state(:,5);
azi  = solution.phase.state(:,6);

aoa  = solution.phase.state(:,7);
bank  = solution.phase.state(:,8);
mFuel  = solution.phase.state(:,9);
throttle  = solution.phase.state(:,10);


lonlist(i) = lon0;
latlist(i) = lat0;
jlist(i) = j;

mFuelList(i) = mFuel(1)


i = i+1
%         end
%     end
% end




% Flap_deflection  = solution.phase.control(:,4);

% [raddot,londot,latdot,fpadot,vdot,azidot, q, M, Fd, rho,L,Fueldt,T,trim_constraint] = VehicleModelReturn(fpa, rad, v,auxdata,azi,lat,lon,aoa,bank,throttle, mFuel,Flap_deflection);
[raddot,londot,latdot,fpadot,vdot,azidot, q, M, Fd, rho,L,Fueldt,T,q1] = VehicleModelReturn(fpa, rad, v,auxdata,azi,lat,lon,aoa,bank,throttle, mFuel);


t = solution.phase(1).time;

% throttle  = solution.phase.control(:,3);

m = Stage2.mStruct+mFuel;

forward0 = [alt0,fpa0,speed0,azi0,lat0,lon0, mFuel(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(t,aoa,f_t),ControlInterp(t,bank,f_t),ControlInterp(t,throttle,f_t)),t(1:end),forward0);

altitude  = (solution.phase(1).state(:,1)-auxdata.Re);
figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(t,altitude);

gamma  = solution.phase.state(:,5);
figure(213)
hold on
plot(f_t(1:end),f_y(:,2));
plot(t,gamma);

latitude  = solution.phase.state(:,3);
figure(214)
hold on
plot(f_y(:,6),f_y(:,5));
plot(lon,lat);

figure(215)
hold on
plot(f_t(1:end),f_y(:,7));
plot(t,mFuel);

save('Solution.mat','solution')
% ----------- Reusable Launch Vehicle Entry Example ------------%
% This example is taken verbatim from the following reference:  %
% Betts, J. T., Practical Methods for Optimal Control Using     %
% Nonlinear Programming, SIAM Press, Philadelphia, 2009.        %
% --------------------------------------------------------------%
%close all
clear all
clc

% cft2m = 0.3048;
% cft2km = cft2m/1000;
% cslug2kg = 14.5939029;


mstruct = 4910.5 - 132.8 + 179.41; % mass of everything but fuel from dawids work

%% Inputs ============================================
%Take input of aero
aero = importdata('SPARTANaero.txt');

interp.Cl_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,3));
interp.Cd_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,4));


[MList,AOAList] = ndgrid(unique(aero(:,1)),unique(aero(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid = [];
Cd_Grid = [];

for i = 1:numel(MList)
    M_temp = MList(i);
    AoA_temp = AOAList(i);
    
    Cl_temp = interp.Cl_scattered(M_temp,AoA_temp);
    Cd_temp = interp.Cd_scattered(M_temp,AoA_temp);
    
    I = cell(1, ndims(MList)); 
    [I{:}] = ind2sub(size(MList),i);
    
    Cl_Grid(I{(1)},I{(2)}) = Cl_temp;
    Cd_Grid(I{(1)},I{(2)}) = Cd_temp;

end

auxdata.interp.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
auxdata.interp.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');

% Produce Atmosphere Data
auxdata.Atmosphere = dlmread('atmosphere.txt');
%=============================================== 

%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)
auxdata.S    = 62;                    % Vehicle Reference Area (m^2)
auxdata.cl   = [-0.2070 1.6756];               % Parameters for Lift Coefficient
auxdata.cd   = [0.0785 -0.3529 2.0400];        % Parameters for Drag Coefficient


auxdata.b    = [0.07854 -0.061592 0.00621408]; % Parameters for Heat Rate Model
auxdata.H    = 7254.24;                        % Density Scale Height (m)
auxdata.al   = [-0.20704 0.029244];            % Parameters for Heat Rate Model
auxdata.rho0 = 1.225570827014494;              % Sea Level Atmospheric Density (kg/m^3)
auxdata.mu   = 3.986031954093051e14;           % Earth Gravitational Parameter (m^3/s^2) 

auxdata.mass = mstruct;               % Vehicle Mass (kg)

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0     = 0;
alt0   = 35000;   
rad0   = alt0+auxdata.Re;
altf   = 30000;   
radf   = altf+auxdata.Re;
lon0   = 0;
lat0   = 0;
speed0 = +2872.88;
speedf = +762;
fpa0   = 3*pi/180; 
fpaf   = 0*pi/180;
azi0   = +90*pi/180; 
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
aziMin = -360*pi/180; aziMax =  360*pi/180;
aoaMin = 0;  aoaMax = 10*pi/180;
bankMin = -90*pi/180; bankMax =   90*pi/180;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [rad0, lon0, lat0, speed0, fpa0, azi0];
bounds.phase.initialstate.upper = [rad0, lon0, lat0, speed0, fpa0, azi0];
bounds.phase.state.lower = [radMin, lonMin, latMin, speedMin, fpaMin, aziMin];
bounds.phase.state.upper = [radMax, lonMax, latMax, speedMax, fpaMax, aziMax];
bounds.phase.finalstate.lower = [radf, lonMin, latMin, speedMin, fpaf, azif];
bounds.phase.finalstate.upper = [radf, lonMax, latMax, speedMax, fpaf, azif];
bounds.phase.control.lower = [aoaMin, bankMin];
bounds.phase.control.upper = [aoaMax, bankMax];
% bounds.phase.control.lower = [aoaMin];
% bounds.phase.control.upper = [aoaMax];
%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess              = [0; 500];
radGuess            = [rad0; radf];
lonGuess            = [lon0; lon0+10*pi/180];
latGuess            = [lat0; lat0+10*pi/180];
speedGuess          = [speed0; speedf];
fpaGuess            = [fpa0; fpaf];
aziGuess            = [azi0; azif];
aoaGuess            = [8*pi/180; 8*pi/180];
bankGuess           = [-25*pi/180; -25*pi/180];
guess.phase.state   = [radGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess];
guess.phase.control = [aoaGuess, bankGuess];
% guess.phase.control = [aoaGuess];
guess.phase.time    = tGuess;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 2;
mesh.colpointsmin = 3;
mesh.colpointsmax = 30;
mesh.tolerance    = 1e-6;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @rlvEntryContinuous;
setup.functions.endpoint             = @rlvEntryEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);

solution = output.result.solution;
aoa       = solution.phase(1).control(:,1);
bank      = solution.phase(1).control(:,2);
t = solution.phase(1).time;

m = mstruct;
forward0 = [alt0,fpa0,speed0,azi0-deg2rad(90),lat0,lon0];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata.interp, auxdata.Atmosphere,ControlInterp(t,aoa,f_t),ControlInterp(t,bank,f_t)),t(1:end),forward0);

altitude  = (solution.phase(1).state(:,1)-auxdata.Re);
figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(t,altitude);
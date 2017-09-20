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
%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)
% auxdata.S    = 249.9091776;                    % Vehicle Reference Area (m^2)
% auxdata.cl   = [-0.2070 1.6756];               % Parameters for Lift Coefficient
% auxdata.cd   = [0.0785 -0.3529 2.0400];        % Parameters for Drag Coefficient
% auxdata.b    = [0.07854 -0.061592 0.00621408]; % Parameters for Heat Rate Model
% auxdata.H    = 7254.24;                        % Density Scale Height (m)
% auxdata.al   = [-0.20704 0.029244];            % Parameters for Heat Rate Model
% auxdata.rho0 = 1.225570827014494;              % Sea Level Atmospheric Density (kg/m^3)
% auxdata.mu   = 3.986031954093051e14;           % Earth Gravitational Parameter (m^3/s^2) 
% auxdata.mass = 92079.2525560557;               % Vehicle Mass (kg)

auxdata.Atmosphere = dlmread('atmosphere.txt');
auxdata.Aero = dlmread('AeroCoeffs.txt');
Aero = auxdata.Aero;
auxdata.Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));
% 
auxdata.Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));

auxdata.CP_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,7));

auxdata.CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

auxdata.Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

auxdata.ThirdStagem = 3300;

auxdata.phi0 = -0.13;
auxdata.xi0 = 0;
auxdata.zeta0 = 1.78;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0     = 0;
alt0   = 35000;   
% rad0   = alt0+auxdata.Re;
altf   = 84000;   
radf   = altf+auxdata.Re;
lon0   = 0;
lat0   = 0;
speed0 = 2900;
speedf = 7000;
fpa0   = deg2rad(3); 
fpaf   = 0;
azi0   = +90*pi/180; 
azif   = -90*pi/180;

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
tfMin = 0;            tfMax = 3000;
altMin = 30000;  altMax = 84000;
lonMin = -pi;         lonMax = -lonMin;
latMin = -70*pi/180;  latMax = -latMin;
speedMin = 10;        speedMax = 8000;
fpaMin =deg2rad(-5);  fpaMax =  deg2rad(30);
aziMin = -180*pi/180; aziMax =  180*pi/180;
aoaMin = 0;  aoaMax = deg2rad(20);

aoadotMin = -deg2rad(1);
aoadotMax = deg2rad(1);
% bankMin = -90*pi/180; bankMax =   1*pi/180;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;

% bounds.phase.initialstate.lower = [rad0, lon0, lat0, speed0, fpa0, azi0, 3300];
% bounds.phase.initialstate.upper = [rad0, lon0, lat0, speed0, fpa0, azi0, 3300];
% 
% bounds.phase.state.lower = [radMin, lonMin, latMin, speedMin, fpaMin, aziMin, 0];
% bounds.phase.state.upper = [radMax, lonMax, latMax, speedMax, fpaMax, aziMax, 3300];
% 
% bounds.phase.finalstate.lower = [radf, lonMin, latMin, speedMin, fpaMin, aziMin, 0];
% bounds.phase.finalstate.upper = [radf, lonMax, latMax, speedMax, fpaMax, aziMax, 3300];

bounds.phase.initialstate.lower = [alt0, speed0, fpa0, auxdata.ThirdStagem, aoaMin];
bounds.phase.initialstate.upper = [alt0, speed0, fpa0, auxdata.ThirdStagem, aoaMax];

bounds.phase.state.lower = [altMin,speedMin, fpaMin, 0, aoaMin];
bounds.phase.state.upper = [altMax, speedMax, fpaMax, auxdata.ThirdStagem, aoaMax];

bounds.phase.finalstate.lower = [altMin, speedMin, 0, 0, aoaMin];
bounds.phase.finalstate.upper = [altMax, speedMax, fpaMax, auxdata.ThirdStagem, aoaMax];

bounds.phase.control.lower = [aoadotMin];
bounds.phase.control.upper = [aoadotMax];

bounds.phase.path.lower = [-deg2rad(8), -inf];
bounds.phase.path.upper = [deg2rad(8), 0];

bounds.eventgroup.lower = 100000;
bounds.eventgroup.upper = 566000;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess              = [0; 250];
altGuess            = [alt0; 80000];
lonGuess            = [lon0; lon0];
latGuess            = [lat0; lon0];
speedGuess          = [speed0; speedf];
fpaGuess            = [fpa0; rad2deg(10)];
aziGuess            = [azi0; azif];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(5); deg2rad(20)];
% bankGuess           = [0; 0];
% guess.phase.state   = [radGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, mGuess];
guess.phase.state   = [altGuess, speedGuess, fpaGuess, mGuess, aoaGuess];
guess.phase.control = [0;0];
guess.phase.time    = tGuess;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 2;
mesh.colpointsmin = 3;
mesh.colpointsmax = 20;
mesh.tolerance    = 1e-4;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @ThirdStageContinuous;
setup.functions.endpoint             = @rlvEntryEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.maxiterations = 100;
setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%


output = gpops2(setup);


%%
solution = output.result.solution;

alt  = solution.phase.state(:,1);
v    = solution.phase.state(:,2);
gamma  = solution.phase.state(:,3);
m    = solution.phase.state(:,4);
aoa    = solution.phase.state(:,5);
aoadot       = solution.phase(1).control(:,1);

forward0 = [alt(1),v(1),gamma(1),m(1)];

time = solution.phase(1).time;

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,aoa,f_t)),time(1:end),forward0);

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T, L, D] = ThirdStageDyn(alt,0,0,gamma,v,deg2rad(97),m,aoa,auxdata);

[AltF_actual, vF, Altexo, vexo, timeexo, mpayload, Alpha, mexo,qexo,gammaexo,Dexo,zetaexo,phiexo, incexo,Texo,CLexo,Lexo,inc_diff] = ThirdStageSim(alt(end),gamma(end),v(end), 0,0, deg2rad(97), m(end), auxdata);

figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

figure(213)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,v);


figure(214)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time,m);


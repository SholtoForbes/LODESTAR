%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scramjet Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the DIDO proprietary optimisation software
% startup.m must be run before this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\')
addpath('..\..\')
addpath('..\EngineData')

%% Atmosphere Data %%======================================================
Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = Atmosphere;


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

%Take input of aero
% aero = importdata('SPARTANaero.txt');
aero = importdata('communicator.txt');
interp.Cl_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,3));
interp.Cd_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,4));
% interp.Cm_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,5));

[MList,AOAList] = ndgrid(unique(aero(:,1)),unique(aero(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid = [];
Cd_Grid = [];
% Cm_Grid = [];

for i = 1:numel(MList)
    M_temp = MList(i);
    AoA_temp = AOAList(i);
    
    Cl_temp = interp.Cl_scattered(M_temp,AoA_temp);
    Cd_temp = interp.Cd_scattered(M_temp,AoA_temp);
%     Cm_temp = interp.Cm_scattered(M_temp,AoA_temp);
    
    I = cell(1, ndims(MList)); 
    [I{:}] = ind2sub(size(MList),i);
    
    Cl_Grid(I{(1)},I{(2)}) = Cl_temp;
    Cd_Grid(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid(I{(1)},I{(2)}) = Cm_temp;

end

auxdata.interp.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
auxdata.interp.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
% auxdata.interp.Cm_spline = griddedInterpolant(MList,AOAList,Cm_Grid,'spline','nearest');

%Take input of aero
% flapaero = importdata('SPARTAN_Flaps.txt');
% 
% auxdata.interp.flap_momentCl_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,3), 'linear', 'nearest');
% auxdata.interp.flap_momentCd_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,4), 'linear', 'nearest');
% auxdata.interp.flap_momentdef_scattered = scatteredInterpolant(aero(:,1),aero(:,5),aero(:,2), 'linear', 'nearest');
% 


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


% Primal Bounds
bounds.phase.state.lower = [Stage2.Bounds.Alt(1), Stage2.Bounds.v(1), Stage2.Bounds.gamma(1), Stage2.Bounds.mFuel(1), Stage2.Bounds.zeta(1),0];
bounds.phase.state.upper = [Stage2.Bounds.Alt(2), Stage2.Bounds.v(2), Stage2.Bounds.gamma(2), Stage2.Bounds.mFuel(2), Stage2.Bounds.zeta(2), deg2rad(10)];

% Initial States
bounds.phase.initialstate.lower = [Stage2.Bounds.Alt(1), Stage2.Initial.v, Stage2.Bounds.gamma(1), Stage2.Initial.mFuel, Stage2.Bounds.zeta(1),0] ;
bounds.phase.initialstate.upper = [Stage2.Bounds.Alt(2), Stage2.Initial.v, Stage2.Bounds.gamma(2), Stage2.Initial.mFuel, Stage2.Bounds.zeta(2), deg2rad(10)];

% End States
bounds.phase.finalstate.lower = [Stage2.Bounds.Alt(1), Stage2.Bounds.v(1), Stage2.End.gammaOpt(1), Stage2.End.mFuel, Stage2.End.Zeta,0];
bounds.phase.finalstate.upper = [Stage2.Bounds.Alt(2), Stage2.Bounds.v(2), Stage2.End.gammaOpt(2), Stage2.Initial.mFuel, Stage2.End.Zeta, deg2rad(10)];

% Control Bounds
bounds.phase.control.lower = 0.01;
bounds.phase.control.upper = 0.01; 

% Time Bounds

bounds.phase.initialtime.lower = 0;
bounds.phase.initialtime.upper = 0;
bounds.phase.finaltime.lower = 30;
bounds.phase.finaltime.upper = Stage2.Bounds.time(2);

%% Define Path Constraints
% This limits the dynamic pressure.
if const == 1 || const == 14 || const == 15
    bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [50000];
elseif const == 12
    bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [55000];
elseif const == 13
    bounds.phase.path.lower = [0 ];
    bounds.phase.path.upper = [45000];
elseif const ==3 || const == 32
        bounds.phase.path.lower = [0];
    bounds.phase.path.upper = [50010];
end


%%  Guess =================================================================

guess.phase.state(:,1)   = [24000;24000];
guess.phase.state(:,2)   = Stage2.Guess.v.';
guess.phase.state(:,3)   = Stage2.Guess.gamma.';
guess.phase.state(:,4) 	= [Stage2.Guess.mFuel(1);Stage2.Guess.mFuel(1)];
guess.phase.state(:,5)   = Stage2.Guess.zeta.';
guess.phase.state(:,6)   = [deg2rad(5);deg2rad(5)];
guess.phase.control      = Stage2.Guess.control.';
guess.phase.time          = Stage2.Guess.time.';


%% Start Iterative Plot
% figure(10)
% plot(linspace(guess.time(1),guess.time(2),algorithm.nodes),linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes),'k')
% iterative_V(end+1,:) = linspace(guess.states(1,1),guess.states(1,2),algorithm.nodes);
% iterative_t(end+1,:) = linspace(guess.time(1),guess.time(2),algorithm.nodes);
% 
% filename = 'testnew51.gif';
% frame = getframe(10);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif', 'Loopcount',0);


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 1;
mesh.colpointsmin = 10;
mesh.colpointsmax = 30;
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
setup.nlp.ipoptoptions.maxiterations = 500;
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%


output = gpops2(setup);

%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================
% Assign the primal variables
alt = output.result.solution.phase.state(:,1).';
v = output.result.solution.phase.state(:,2).'; 
gamma = output.result.solution.phase.state(:,3).'; 
mfuel = output.result.solution.phase.state(:,4).'; 
% gammadot = output.result.solution.phase.state(:,5).';
zeta = output.result.solution.phase.state(:,5).';
alpha = output.result.solution.phase.state(:,6).';

omegadot  = output.result.solution.phase.control.'; 

time = output.result.solution.phase.time.';

[altdot, dfuel, Fueldt, vdot, q, M, D, Thrust, flapdeflection, rho,zeta,phi,eq,zetadot,xi, gammadot] = VehicleModel(time, gamma, alt, v, mfuel,interp,const, interp.Atmosphere,zeta,auxdata,alpha,gamma);


% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point



% cd('../ThirdStage')
% [ThirdStagePayloadMass,ThirdStageControls,ThirdStageZeta,ThirdStagePhi,ThirdStageAlt,ThirdStagev,ThirdStaget,ThirdStageAlpha,ThirdStagem,ThirdStagegamma,ThirdStageq] = ThirdStageOptm(V(end),gamma(end),v(end), phi(end),zeta(end), 1);
% ThirdStagePayloadMass
% cd('../SecondStage')
ThirdStagePayloadMass = 0;



%% Check End Point
if v(end) < min(ThirdStageData(:,5)) || v(end) > max(ThirdStageData(:,5)) || gamma(end) < min(ThirdStageData(:,4)) || gamma(end) > max(ThirdStageData(:,4)) || alt(end) < min(ThirdStageData(:,3)) || alt(end) > max(ThirdStageData(:,3))
msgbox('ERROR: OPTIMISATION INVALID. End point outside of third stage matrix.');
end

if abs(output.result.objective - ThirdStagePayloadMass) > 1;
msgbox('ERROR: Check end point. Third stage payload mass calculated from optimised trajectory is significantly different from interpolated result.');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt)


dt = time(2:end)-time(1:end-1); % Time change between each node pt
FuelUsed = zeros(1,nodes-1);
FuelUsed(1) = dt(1)*Fueldt(1);
for i = 2:nodes-1
    FuelUsed(i) = dt(i).*Fueldt(i) + FuelUsed(i-1);
end


% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
end

% Separation_LD = lift(end)/Fd(end)
% 
figure(201)

subplot(5,5,[1,10])
hold on
plot(H, alt)
% plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')

dim = [.7 .52 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1000 - mfuel(end)) ' kg']},'FitBoxToText','on');  


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
plot(time, mfuel + 8755.1 - 994)
title('Vehicle Mass (kg)')



subplot(5,5,17)
plot(time, Thrust)
title('Thrust (N)')

Isp = Thrust./Fueldt./9.81;
IspNet = (Thrust-D)./Fueldt./9.81;

subplot(5,5,18)
plot(time, Isp)
title('Isp')

subplot(5,5,19)
plot(time, IspNet)
title('Net Isp')

% subplot(5,5,20)
% plot(time, flapdeflection)
% title('Flap Deflection (deg)')

subplot(5,5,21)
plot(time, rad2deg(alpha))
title('Angle of Attack (deg)')

% subplot(5,5,22);
% plot(time, dual.dynamics);
% title('costates')
% xlabel('time');
% ylabel('costates');
% legend('\lambda_1', '\lambda_2', '\lambda_3');
% 
% subplot(5,5,23)
% Hamiltonian = dual.Hamiltonian(1,:);
% plot(time,Hamiltonian);
% title('Hamiltonian')

subplot(5,5,24)
hold on
plot(time, rad2deg(gammadot))
title('Trajectory Angle Change Rate (Deg/s)')

subplot(5,5,25)
hold on
plot(time, rad2deg(omegadot))
title('Omegadot Control (Deg/s2)')


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
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mfuel(1) - mfuel(end)) ' kg']},'FitBoxToText','on');  

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
line(time, [Alpha(1:end-1) Alpha(end-1)],'Parent',ax3,'Color','k', 'LineStyle','-')

line(time, [flapdeflection(1:end-1) flapdeflection(end-1)],'Parent',ax3,'Color','k', 'LineStyle','--')


% line(time, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
line(time, eq.*10,'Parent',ax3,'Color','k', 'LineStyle','-.')

line(time, IspNet./(10^2),'Parent',ax3,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax3, 'AoA (degrees)','Flap Deflection (degrees)', 'Equivalence Ratio x 10', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

saveas(figure(202),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'SecondStage.fig']);

% 
% dat_temp1 = get(ax1,'children');
% fig_temp1 = figure;
% ax_temp1 = axes;
% temp_fig1 = copyobj(dat_temp1,ax_temp1);
% title('Trajectory')
% xlabel('Earth Normal Distance Flown (km)')
% ylabel('Vertical Position (km)')
% dim = [.55 .15 .2 .2];
% annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mfuel(1) - mfuel(end)) ' kg']},'FitBoxToText','on');  
%  set(fig_temp1, 'Position', [100, 100, 900, 400]);
%  saveas(fig_temp1,[pwd sprintf('../ArchivedResults/%s/FlightPath',Timestamp)]);
%  close fig_temp1;
%  
% dat_temp2 = get(ax2,'children');
% fig_temp2 = figure;
% ax_temp2 = axes;
% temp_fig2 = copyobj(dat_temp2,ax_temp2);
%  set(fig_temp2, 'Position', [100, 100, 900, 400]);
%  h = legend(ax_temp2,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
% rect1 = [0.22, 0.85, .25, .25];
% % set(h, 'Position', rect1)
% set(h,'location','bestoutside')
% xlabel('time (s)')
 
 
% temp_fig3 = copyobj(sp3,ax3);
%  set(FigHandle, 'Position', [100, 100, 900, 400]);

figure(203)

subplot(2,5,[1,5]);

line(time, dual.dynamics(1,:),'Color','k', 'LineStyle','-');
line(time, dual.dynamics(2,:),'Color','k', 'LineStyle','--');
line(time, dual.dynamics(3,:),'Color','k', 'LineStyle','-.');
line(time, dual.dynamics(4,:),'Color','k', 'LineStyle',':');
line(time, dual.dynamics(5,:),'Color','k', 'LineStyle','-','LineWidth',2);
title('costates')
xlabel('time');
ylabel('Costates');
% axis([0,time(end),-1,1])
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4','\lambda_5');

subplot(2,5,[6,10])
Hamiltonian = dual.Hamiltonian(1,:);
plot(time,Hamiltonian,'Color','k');
axis([0,time(end),-1,1])
title('Hamiltonian')


% save results
dlmwrite('primal.txt', [primal.states;primal.controls;primal.nodes;q;IspNet;Alpha;M;eq;flapdeflection;phi]);
dlmwrite('payload.txt', ThirdStagePayloadMass);
dlmwrite('dual.txt', [dual.dynamics;dual.Hamiltonian]);
dlmwrite('ThirdStage.txt',[ThirdStageZeta;ThirdStagePhi;ThirdStageAlt;ThirdStagev;ThirdStaget;[ThirdStageAlpha 0];ThirdStagem;ThirdStagegamma;[ThirdStageq 0]]);
dlmwrite('LD.txt', Separation_LD);


copyfile('primal.txt',sprintf('../ArchivedResults/%s/primal_%s.txt',Timestamp,Timestamp))
copyfile('dual.txt',sprintf('../ArchivedResults/%s/dual_%s.txt',Timestamp,Timestamp))
copyfile('payload.txt',sprintf('../ArchivedResults/%s/payload_%s.txt',Timestamp,Timestamp))
copyfile('ThirdStage.txt',sprintf('../ArchivedResults/%s/ThirdStage_%s.txt',Timestamp,Timestamp))
copyfile('LD.txt',sprintf('../ArchivedResults/%s/LD_%s.txt',Timestamp,Timestamp))
primal_old = primal;

ts = timeseries(Isp,time);
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

plot(time,dLHdu,time,mu_1,time,mu_2,time,mu_3,time,mu_4,time,mu_5,time,mu_u);
legend('dLHdu','mu_1','mu_2','mu_3','mu_4','mu_5','mu_u');
title('Validation')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD INTEGRATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This simply tests that the system dynamics hold, as the
% Pseudospectral method may not converge to a realistic
% solution



gamma_F = cumtrapz(time,gammadot)+ gamma(1);

gammadot_F = cumtrapz(time,omegadot) + gammadot(1);


v_F = cumtrapz(time,a);
v_F = v_F + v(1);

V_F = cumtrapz(time,v_F.*sin(gamma_F));
V_F = V_F + alt(1);

mfuel_F = cumtrapz(time,-Fueldt);
mfuel_F = mfuel_F + mfuel(1);

figure(206)

subplot(5,1,1)
plot(time,gamma_F,time,gamma);
title('Forward Simulation Comparison');
% note this is just a trapezoidal rule check, may not be exactly accurate
subplot(5,1,2)
plot(time,v_F,time,v);
subplot(5,1,3)
plot(time,V_F,time,alt);
subplot(5,1,4)
plot(time,mfuel_F,time,mfuel);
subplot(5,1,4)
plot(time,gammadot_F,time,gammadot);

% Compute difference with CADAC for constant dynamic pressure path

t_diff = time - [0 time(1:end-1)];
if const == 3
    CADAC_DATA = dlmread('TRAJ.ASC');
    CADAC_Alpha = interp1(CADAC_DATA(:,2),CADAC_DATA(:,4),M(1:68)); % 1:68 gives mach numbers that align, may need to change this
    CADAC_V = interp1(CADAC_DATA(:,2),CADAC_DATA(:,11),M(1:68));
    MeanError_V = sum(abs(CADAC_V - alt(1:68))./alt(1:68).*t_diff(1:68))/time(end)
    MeanError_Alpha = sum(abs(CADAC_Alpha - Alpha(1:68))./Alpha(1:68).*t_diff(1:68))/time(end)
end

% if PayloadGrid(phi(end),zeta(end),V(end)+10,gamma(end),v(end)) - PayloadGrid(phi(end),zeta(end),V(end),gamma(end),v(end)) < 0
%     disp('Check Third Stage Payload Matrix, May Have Found False Maxima')
% end
if PayloadGrid(alt(end)+10,gamma(end),v(end)) - PayloadGrid(alt(end),gamma(end),v(end)) < 0
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

forward0 = [alt(1),phi(1),gamma(1),v(1),zeta(1),Stage2.mStruct+Stage3.mTot+Stage2.mFuel];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(time,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,interp),time,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(time,alpha,f_t),communicator,communicator_trim,interp.Atmosphere,const,auxdata.interp,AlphaInterp(time,lift,f_t),AlphaInterp(time,Fd,f_t),AlphaInterp(time,Thrust,f_t),AlphaInterp(time,flapdeflection,f_t)),time(1:end),forward0);

figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,phi);

subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time,gamma);

subplot(7,1,5)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time,v);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,5));
plot(time,zeta);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,6));
plot(time,Stage2.mStruct+Stage3.mTot+mfuel);

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






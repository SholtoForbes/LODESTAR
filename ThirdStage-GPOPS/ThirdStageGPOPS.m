% function mpayload = ThirdStageGPOPS(alt0,gamma0,v0, phi0, zeta0);

alt0 = 35000;
gamma0 = deg2rad(3);
v0 = 2850;
phi0 = -0.13;
zeta0 = 1.78


%-------------------------------------------------------------------------%
%------------------ Provide Auxiliary Data for Problem -------------------%
%-------------------------------------------------------------------------%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)

interp.Atmosphere = dlmread('atmosphere.txt');

auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data

auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data

auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data

auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 

auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 

auxdata.Aero = dlmread('AeroCoeffs.txt');
Aero = auxdata.Aero;
auxdata.Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));
% 
auxdata.Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));

auxdata.CP_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,7));

auxdata.CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

auxdata.Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

auxdata.ThirdStagem = 3300;

auxdata.phi0 = phi0;
auxdata.xi0 = 0;
auxdata.zeta0 = zeta0;
%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0     = 0;
% altf   = 84000;   
% radf   = altf+auxdata.Re;
% v0 = 2900;
% vf = 7000;
% gammaf   = 0;

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
tfMin = 7;            tfMax = 200;
altMin = alt0;  altMax = 100000;
phiMin = -0.2;         phiMax = 0.1;
vMin = 10;        vMax = 8000;
gammaMin =deg2rad(-5);  gammaMax =  deg2rad(30);
zetaMin = deg2rad(90); zetaMax =  deg2rad(110);
aoaMin = 0;  aoaMax = deg2rad(20);

aoadotMin = -deg2rad(1);
aoadotMax = deg2rad(1);

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;


bounds.phase.initialstate.lower = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMin, phi0, zeta0];
bounds.phase.initialstate.upper = [alt0, v0, gamma0, auxdata.ThirdStagem, aoaMax, phi0, zeta0];

bounds.phase.state.lower = [altMin,vMin, gammaMin, 0, aoaMin, phiMin, zetaMin];
bounds.phase.state.upper = [altMax, vMax, gammaMax, auxdata.ThirdStagem, aoaMax, phiMax, zetaMax];

bounds.phase.finalstate.lower = [altMin, vMin, 0, 0, 0, phiMin, zetaMin];
bounds.phase.finalstate.upper = [altMax, vMax, gammaMax, auxdata.ThirdStagem, 0, phiMax, zetaMax];

bounds.phase.control.lower = [aoadotMin];
bounds.phase.control.upper = [aoadotMax];

bounds.phase.path.lower = [-deg2rad(8), -inf];
bounds.phase.path.upper = [deg2rad(8), 0];

bounds.eventgroup.lower = 100000;
bounds.eventgroup.upper = 566000;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess              = [0; 150];
altGuess            = [alt0; 90000];
vGuess          = [v0; 7000];
gammaGuess            = [gamma0; rad2deg(10)];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(20); deg2rad(20)];
phiGuess = [phi0;-0.05];
zetaGuess = [zeta0;zeta0];
guess.phase.state   = [altGuess, vGuess, gammaGuess, mGuess, aoaGuess, phiGuess, zetaGuess];
guess.phase.control = [0;0];
guess.phase.time    = tGuess;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 4
mesh.colpointsmax = 30;
mesh.tolerance    = 1e-5;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @ThirdStageContinuous;
setup.functions.endpoint             = @ThirdStageEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 0;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.maxiterations = 500;
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
phi    = solution.phase.state(:,6);
zeta    = solution.phase.state(:,7);
aoadot       = solution.phase(1).control(:,1);

forward0 = [alt(1),v(1),gamma(1),m(1),phi(1),zeta(1)];

time = solution.phase(1).time;

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,aoa,f_t),ControlInterp(time,aoadot,f_t)),time(1:end),forward0);

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T, L, D, q] = ThirdStageDyn(alt,gamma,v,m,aoa,time,auxdata,aoadot,phi,zeta);


xi(1) = auxdata.xi0;
for i = 2:length(time)
    xi(i) = xi(i-1) + xidot(i-1)*(time(i)-time(i-1));
end

[AltF_actual, vF, altexo, vexo, timeexo, mpayload, Alpha, mexo,qexo,gammaexo,Dexo,zetaexo,phiexo, incexo,Texo,CLexo,Lexo,inc_diff] = ThirdStageSim(alt(end),gamma(end),v(end), phi(end),xi(end), zeta(end), m(end), auxdata);



figure(212)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

figure(213)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,v);


%% Plotting
% if plotflag == 1
figure(214)
    addpath('addaxis')
    hold on

    plot([time; timeexo.'+time(end)], [alt; altexo.']/1000, 'LineStyle', '-','Color','k', 'lineWidth', 2.2)
    plot([time; timeexo.'+time(end)],[q;qexo.';qexo(end)]/1000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
    plot([time; timeexo.'+time(end)],[rad2deg(aoa);0*ones(length(timeexo),1)], 'LineStyle', '--','Color','k', 'lineWidth', 0.7)
    ylabel('Altitude (km), Dynamic Pressure (kPa), Angle of Attack (deg)');
    

    addaxis([time; timeexo.'+time(end)],[v;vexo.'], [0 7000], 'LineStyle', '--','Color','k', 'lineWidth', 1.8)
    addaxisplot([time; timeexo.'+time(end)],[ m;mexo.';mexo(end)],2, 'LineStyle', ':','Color','k', 'lineWidth', 1.3)
    addaxislabel(2,'Velocity (m/s), Mass (kg)');


    addaxis([time; timeexo.'+time(end)],[rad2deg(Vec_angle);0*ones(length(timeexo),1)], 'LineStyle', ':','Color','k', 'lineWidth', 2.1)
    addaxisplot([time; timeexo.'+time(end)], [rad2deg(gamma);rad2deg(gammaexo).'],3, 'LineStyle', '-','Color','k', 'lineWidth', .6)
    addaxislabel(3,'Thrust Vector Angle (deg), Trajectory Angle (deg)');

    legend(  'Altitude','Dynamic Pressure','Angle of Attack', 'Velocity',  'Mass', 'Thrust Vector Angle', 'Trajectory Angle' );
    xlabel('Time (s)');
    xlim([0 timeexo(end)+time(end)])
    box off
    % Write data to file
    dlmwrite('ThirdStageData',[time, alt, v, m,q ,gamma,D ,zeta], ' ')

%     Integrated_Drag = cumtrapz(time(1:end-1),D) ;
%     Integrated_Drag(end)
% end
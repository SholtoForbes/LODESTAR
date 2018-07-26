%%
% Created by Sholto Forbes 25/7/18
% This file processes the GPOPS-2 solution output for the
% rocket-scramjet-rocket launch trajectory.
%%

function [] = Plotter(output,auxdata,mode,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,mRocket,mSpartan,mFuel,h0,v0,bounds)

Timestamp = datestr(now,30)
mkdir('../ArchivedResults', strcat(Timestamp, 'mode', num2str(mode)))


copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',strcat(Timestamp,'mode',num2str(mode))))

save output
movefile('output.mat',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

% Run multiple times if output has mutltiple cells
for j = 1:length(output)




% =========================================================================
% Assign the primal variables
alt21 = output{j}.result.solution.phase(2).state(:,1);
alt22 = output{j}.result.solution.phase(3).state(:,1);
lon21 = output{j}.result.solution.phase(2).state(:,2);
lon22 = output{j}.result.solution.phase(3).state(:,2);
lat21 = output{j}.result.solution.phase(2).state(:,3);
lat22 = output{j}.result.solution.phase(3).state(:,3);
v21 = output{j}.result.solution.phase(2).state(:,4); 
v22 = output{j}.result.solution.phase(3).state(:,4); 
gamma21 = output{j}.result.solution.phase(2).state(:,5); 
gamma22 = output{j}.result.solution.phase(3).state(:,5); 
zeta21 = output{j}.result.solution.phase(2).state(:,6);
zeta22 = output{j}.result.solution.phase(3).state(:,6);
alpha21 = output{j}.result.solution.phase(2).state(:,7);
alpha22 = output{j}.result.solution.phase(3).state(:,7);
eta21 = output{j}.result.solution.phase(2).state(:,8);
eta22 = output{j}.result.solution.phase(3).state(:,8);
mFuel21 = output{j}.result.solution.phase(2).state(:,9); 
mFuel22 = output{j}.result.solution.phase(3).state(:,9); 

throttle22 = output{j}.result.solution.phase(3).state(:,10);

aoadot21  = output{j}.result.solution.phase(2).control(:,1); 
etadot21  = output{j}.result.solution.phase(2).control(:,2); 

aoadot22  = output{j}.result.solution.phase(3).control(:,1); 
etadot22  = output{j}.result.solution.phase(3).control(:,2); 

time21 = output{j}.result.solution.phase(2).time;
time22 = output{j}.result.solution.phase(3).time;

figure(01)
subplot(9,1,1)
hold on
plot(time21,alt21)
plot(time22,alt22)
subplot(9,1,2)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,3)
hold on
plot(time21,lon21)
plot(time22,lon22)
subplot(9,1,4)
hold on
plot(time21,lat21)
plot(time22,lat22)
subplot(9,1,5)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,6)
hold on
plot(time21,gamma21)
plot(time22,gamma22)
subplot(9,1,7)
hold on
plot(time21,ones(1,length(time21)))
plot(time22,throttle22)


figure(230)
hold on
plot3(lon21,lat21,alt21)
plot3(lon22,lat22,alt22)

  
%     figure(2301)
% hold on
% 
% axesm('pcarree','Origin',[0 rad2deg(lon21(1)) 0])
% geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
% % plotm(rad2deg(lat),rad2deg(lon+lon0))
% plotm(rad2deg(lat21),rad2deg(lon21),'b')
% plotm(rad2deg(lat22),rad2deg(lon22),'r')
%     
%     cities = shaperead('worldcities', 'UseGeoCoords', true);
% lats = extractfield(cities,'Lat');
% lons = extractfield(cities,'Lon');
% geoshow(lats, lons,...
%         'DisplayType', 'point',...
%         'Marker', 'o',...
%         'MarkerEdgeColor', 'r',...
%         'MarkerFaceColor', 'r',...
%         'MarkerSize', 2)

% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point

ThirdStagePayloadMass = -output{j}.result.objective;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt21)

[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho21,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,1, mFuel21,mFuel21(1),mFuel21(end), 1, 0);
[~,~,~,~,~,~, q22, M22, Fd22, rho22,L22,Fueldt22,T22,Isp22,q122,flapdeflection22,heating_rate22] = VehicleModelCombined(gamma22, alt22, v22,auxdata,zeta22,lat22,lon22,alpha22,eta22,throttle22, mFuel22,0,0, 0, 0);

throttle22(M22<5.0) = 0; % remove nonsense throttle points
throttle22(q122<20000) = 0; % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
    
Isp22(M22<5.0) = 0; % remove nonsense throttle points

% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v21(i)*(time21(i+1) - time21(i))*cos(gamma21(i)) + H(i);
end

% Separation_LD = lift(end)/Fd(end)

%% plot Ascent
figure(211)
fig = gcf;
set(fig,'Position',[200 0 850 1200])

subplot(7,2,1)
hold on
plot(time21-time21(1), alt21/1000,'Color','k')
title('Trajectory (km)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

dim = [.55 .7 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Fuel Used: ' num2str(1562 - mFuel21(end)) ' kg']},'FitBoxToText','on');  

subplot(7,2,3)
hold on
plot(time21-time21(1), v21,'Color','k')
title('Velocity (m/s)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,4)
plot(time21-time21(1), M21,'Color','k')
title('Mach no','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,5)
plot(time21-time21(1), q21/1000,'Color','k')
title('Dynamic Pressure (kpa)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,6)
hold on
plot(time21-time21(1), rad2deg(gamma21),'Color','k')
title('Trajectory Angle (Deg)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,7)
plot(time21-time21(1), rad2deg(alpha21),'Color','k')
title('Angle of Attack (deg)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,8)
hold on
plot(time21-time21(1), rad2deg(eta21),'Color','k')
title('Bank Angle (deg)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,9)
plot(time21-time21(1), flapdeflection21,'Color','k')
title('Flap Deflection (deg)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

% Isp1 = T1./Fueldt1./9.81;
IspNet21 = (T21-Fd21)./Fueldt21./9.81;

subplot(7,2,13)
plot(time21-time21(1), T21/1000,'Color','k')
title('Thrust (kN)','FontSize',9)
set(gca,'xticklabels',[])
xlabel('Time (s)','FontSize',9);
xlim([0 time21(end)-time21(1)]);

subplot(7,2,12)
plot(time21-time21(1), IspNet21,'Color','k')
title('Net Isp (s)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,11)
plot(time21-time21(1), mFuel21,'Color','k')
title('Fuel Mass (kg)','FontSize',9)
% xlabel('Time (s)','FontSize',9);
xlim([0 time21(end)-time21(1)]);
set(gca,'xticklabels',[])

subplot(7,2,10)
plot(time21-time21(1), rad2deg(zeta21),'Color','k')
title('Heading Angle (deg)','FontSize',9)
set(gca,'xticklabels',[])
xlim([0 time21(end)-time21(1)]);

subplot(7,2,14)
plot(time21-time21(1), L21./Fd21,'Color','k')
title('L/D','FontSize',9)
xlabel('Time (s)','FontSize',9);
xlim([0 time21(end)-time21(1)]);

ThirdStageFlag = [ones(length(time21),1); zeros(length(time22),1)];

SecondStageStates = [[time21; time22] [alt21; alt22] [lon21; lon22] [lat21; lat22] [v21; v22] [gamma21; gamma22] [zeta21; zeta22] [alpha21; alpha22] [eta21; eta22] [mFuel21; mFuel22] ThirdStageFlag];
dlmwrite(strcat('SecondStageStates',num2str(j)),['time (s) ' 'altitude (m) ' 'longitude (rad) ' 'latitude (rad) ' 'velocity (m/s) ' 'trajectory angle (rad) ' 'heading angle (rad) ' 'angle of attack (rad) ' 'bank angle (rad) ' 'fuel mass (kg) ' 'Third Stage (flag)'],'');
dlmwrite(strcat('SecondStageStates',num2str(j)),SecondStageStates,'-append','delimiter',' ');
movefile(strcat('SecondStageStates',num2str(j)),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

%% Plot Return
figure(221)
fig = gcf;
set(fig,'Position',[200 0 850 1200])

subplot(7,2,1)
xlim([0 time22(end)-time22(1)]);
hold on
plot(time22-time22(1), alt22/1000,'Color','k')
title('Trajectory (km)','FontSize',9)
set(gca,'xticklabels',[])

dim = [.55 .7 .2 .2];
annotation('textbox',dim,'string',{['Fuel Used: ' num2str(mFuel22(1)) ' kg']},'FitBoxToText','on');  

subplot(7,2,3)
xlim([0 time22(end)-time22(1)]);
hold on
plot(time22-time22(1), v22,'Color','k')
title('Velocity (m/s)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,4)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), M22,'Color','k')
title('Mach no','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,5)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), q22/1000,'Color','k')
title('Dynamic Pressure (kpa)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,6)
xlim([0 time22(end)-time22(1)]);
hold on
plot(time22-time22(1), rad2deg(gamma22),'Color','k')
title('Trajectory Angle (Deg)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,7)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), rad2deg(alpha22),'Color','k')
title('Angle of Attack (deg)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,8)
xlim([0 time22(end)-time22(1)]);
ylim([rad2deg(min(eta22))-1 rad2deg(max(eta22))+1])
hold on
plot(time22-time22(1), rad2deg(eta22),'Color','k')
title('Bank Angle (deg)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,9)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), flapdeflection22,'Color','k')
title('Flap Deflection (deg)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,13)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), T22/1000,'Color','k')
title('Thrust (kN)','FontSize',9)

xlabel('Time (s)','FontSize',9);

subplot(7,2,12)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), Isp22,'Color','k')
title('Potential Isp (s)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,11)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), mFuel22,'Color','k')
title('Fuel Mass (kg)','FontSize',9)
set(gca,'xticklabels',[])

subplot(7,2,14)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), throttle22,'Color','k')
title('Throttle','FontSize',9)
xlabel('Time (s)','FontSize',9);

subplot(7,2,10)
hold on
xlim([0 time22(end)-time22(1)]);
plot(time22-time22(1), rad2deg(zeta22),'Color','k')
title('Heading Angle (deg)','FontSize',9)
set(gca,'xticklabels',[])

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
input_test = output{j}.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

H1 = [];
H2 = [];

lambda1 = output{j}.result.solution.phase(2).costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output{j}.result.solution.phase(3).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(3).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end
time21(1:end-1)
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
for i = 1:length(output{j}.result.solution.phase(2).state(1,:))
plot(time21(1:end-1),([diff(output{j}.result.solution.phase(2).state(:,i))./diff(output{j}.result.solution.phase(2).time)] - phaseout_test(2).dynamics(1:end-1,i))./max(abs(phaseout_test(2).dynamics(1:end-1,i))),'--');
end
for i = 1:length(output{j}.result.solution.phase(3).state(1,:))
    if i<= 7 % Plot different line styles when no. of colours exceeded
    plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))));
    else
    plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))),':');
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

alt3  = output{j}.result.solution.phase(4).state(:,1);
v3    = output{j}.result.solution.phase(4).state(:,2);
gamma3  = output{j}.result.solution.phase(4).state(:,3);
m3    = output{j}.result.solution.phase(4).state(:,4);
aoa3    = output{j}.result.solution.phase(4).state(:,5);
lat3    = output{j}.result.solution.phase(4).state(:,6);
zeta3    = output{j}.result.solution.phase(4).state(:,7);
aoadot3       = output{j}.result.solution.phase(4).control(:,1);

forward0 = [alt3(1),v3(1),gamma3(1),m3(1),lat3(1),zeta3(1)];

time3 = output{j}.result.solution.phase(4).time;

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel3_forward(f_t, f_y,auxdata,ControlInterp(time3,aoa3,f_t),ControlInterp(time3,aoadot3,f_t)),time3(1:end),forward0);

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDyn(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);


lon3(1) = lon21(end);
for i = 2:length(time3)
    lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end

[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);




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
    title('Altitude (km','FontSize',9);
    plot([time3; timeexo.'+time3(end)], [alt3; altexo.']/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,2)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[q3;qexo.';qexo(end)]/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,3)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,4)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[v3;v3exo.'],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,5)
    hold on
    title('Mass (kg','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[ m3;mexo.';mexo(end)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,6)
    hold on
    title('Thrust Vector Angle (deg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[rad2deg(Vec_angle3);0*ones(length(timeexo),1)],'Color','k')
    xlabel('Time (s)','FontSize',9);
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,7)
    hold on
    title('Trajectory Angle (deg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)], [rad2deg(gamma3);rad2deg(gammaexo).'],'Color','k')

    xlabel('Time (s)','FontSize',9);
    xlim([time3(1) timeexo(end)+time3(end)])

    
    % Write data to file
    dlmwrite(strcat('ThirdStageStates',num2str(j)),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg) ' 'dynamic pressure (Pa)' 'trajectory angle (rad) ' 'Lift (N)' 'Drag (N)' 'heading angle (rad) ' 'latitude (rad) ' 'angle of attack (rad) '],'');
    dlmwrite(strcat('ThirdStageStates',num2str(j)),[[time3; time3(end)+timeexo'], [alt3; altexo'], [v3; v3exo'], [m3; mexo'; mexo(end)],[q3; qexo'; qexo(end)] ,[gamma3; gammaexo'],[L3; Lexo'; Lexo(end)],[D3; Dexo'; Dexo(end)] ,[zeta3; zetaexo'], [lat3; latexo'], [aoa3; zeros(length(timeexo),1)]],'-append','delimiter',' ')
movefile(strcat('ThirdStageStates',num2str(j)),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

% 
%% First Stage =========================================================

time1 = output{j}.result.solution.phase(1).time.';
% 
alt1 = output{j}.result.solution.phase(1).state(:,1).';
v1 = output{j}.result.solution.phase(1).state(:,2).';
m1 = output{j}.result.solution.phase(1).state(:,3).';
gamma1 = output{j}.result.solution.phase(1).state(:,4).';
alpha1 = output{j}.result.solution.phase(1).state(:,5).';
zeta1 = output{j}.result.solution.phase(1).state(:,6).';
lat1 = output{j}.result.solution.phase(1).state(:,8).';
lon1 = output{j}.result.solution.phase(1).state(:,9).';
% 

% 
FirstStageSMF = (mRocket - mFuel)/(m1(1) - mSpartan);
% 

FirstStageStates = [time1' alt1' v1' m1' gamma1' alpha1' zeta1' lat1' lon1'];

dlmwrite(strcat('FirstStageStates',num2str(j)),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg)' 'trajectory angle (rad) ' 'angle of attack (rad) ' 'heading angle (rad) ' 'latitude (rad)'],'');
dlmwrite(strcat('FirstStageStates',num2str(j)),FirstStageStates,'-append','delimiter',' ');
movefile(strcat('FirstStageStates',num2str(j)),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));


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
tspan = time1; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [h0 v0 m1(1) deg2rad(89.9) lat1(1) zeta1(1)];

[t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t_f,time1,zeta1),ControlFunction(t_f,time1,alpha1),phase,interp,Throttle,Vehicle,Atmosphere), tspan, postpitch0_f);

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
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(y(:,4).') rad2deg(gamma1)],'color','k');
subplot(4,2,2)
hold on
title('Velocity (m/s)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,2).' v1],'color','k');
subplot(4,2,3)
hold on
title('Altitude (km)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,1).'/1000 alt1/1000],'color','k');
subplot(4,2,4)
hold on
title('Angle of Attack (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha1)],'color','k');
subplot(4,2,5)
hold on
title('Mass (kg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,3).' m1],'color','k');
subplot(4,2,6)
hold on
title('Heading Angle (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(y(:,6).') rad2deg(zeta1)],'color','k');
subplot(4,2,7)
hold on
title('Latitude (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(lat1(1)+y(:,8).') rad2deg(lat1)],'color','k');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')
xlabel('Time (s)')
xlim([0,time1(end)+t_prepitch(end)]);

saveas(figure(111),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('FirstStage',num2str(j),'.fig')]);



%% Create Easy Latex Inputs

dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\PayloadToOrbitMode', num2str(mode) ,'}{ ', num2str(round(ThirdStagePayloadMass,1),'%.1f') , '}'), 'delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\12SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(1)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc')

dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\FirstStageSMFMode', num2str(mode) ,'}{ ', num2str(round(FirstStageSMF,3),'%.3f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\23SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(end)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\23SeparationvMode', num2str(mode) ,'}{ ', num2str(round(v21(end),0)) , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\23SeparationqMode', num2str(mode) ,'}{ ', num2str(round(q21(end)/1000,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\23SeparationLDMode', num2str(mode) ,'}{ ', num2str(round(L21(end)/Fd21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\2FlightTimeMode', num2str(mode) ,'}{ ', num2str(round(time21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

qlt20 = find(q3<20000);
dlmwrite(strcat('LatexInputs',num2str(j),'.txt'),strcat('\newcommand{\3qOver20Mode', num2str(mode) ,'}{ ', num2str(round(time3(qlt20(1))-time3(1),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');



%% Bound Check
% Peform check to see if any of the states are hitting their bounds. This
% is an error if the bound is not intended to constrain the state. Fuel
% mass and throttle are not checked, as these will always hit bounds. 

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))-1
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(3).state(1,:))-2
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting upper bound'))
    end
end

% Angle of attack is not checked on third stage, because angle of attack is hard constrained and should be checked manually. 
for i = [1:3 6: length(output{j}.result.solution.phase(4).state(1,:))]
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting upper bound'))
    end
end

%% Plot combined trajectory




        figure(2301)
           fig = gcf;
        set(fig,'Position',[200 0 1000 800])
hold on
axesm('pcarree','Origin',[lat22(1) rad2deg(lon22(1)) 0],'MapLatLimit',[min(rad2deg(lat22))-50 max(rad2deg(lat22))+50],'MapLonLimit',[min(rad2deg(lon22))-50 max(rad2deg(lon21))+50])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
view([-30 35])
% plotm(rad2deg(lat),rad2deg(lon+lon0))
plot3m(rad2deg(lat1),rad2deg(lon1),alt1/1000000,'color',[0 0.8 0],'LineWidth',1.7)
plot3m(rad2deg(lat21),rad2deg(lon21),alt21/1000000,'r','LineWidth',1.7)
plot3m(rad2deg(lat22),rad2deg(lon22),alt22/1000000,'color',[1 0 1],'LineWidth',1.7)
plot3m(rad2deg(lat3),rad2deg(lon3)',alt3/1000000,'b','LineWidth',1.7)
plot3m(rad2deg(latexo),rad2deg(lonexo),altexo/1000000,'c','LineWidth',1.7)

plotm(rad2deg(lat1),rad2deg(lon1),'color',[0 0.8 0])
plotm(rad2deg(lat21),rad2deg(lon21),'r')
plotm(rad2deg(lat22),rad2deg(lon22),'color',[1 0 1])
plotm(rad2deg(lat3),rad2deg(lon3)','b')
plotm(rad2deg(latexo),rad2deg(lonexo),'c')
   
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)

    
for i = 0:20/(time1(end)-time1(1)):1
    time1_temp = (time1(end)-time1(1))*i + time1(1);
    alt1_temp = interp1(time1,alt1,time1_temp);
    lat1_temp = interp1(time1,lat1,time1_temp);
    lon1_temp = interp1(time1,lon1,time1_temp);
    plot3m(rad2deg([lat1_temp  ,lat1_temp ]), rad2deg([lon1_temp,  lon1_temp]),[0,alt1_temp/1000000],'color',[0.3 0.3 0.3]);
end
for i = 0:20/(time21(end)-time21(1)):1
    time21_temp = (time21(end)-time21(1))*i + time21(1);
    alt21_temp = interp1(time21,alt21,time21_temp);
    lat21_temp = interp1(time21,lat21,time21_temp);
    lon21_temp = interp1(time21,lon21,time21_temp);
    plot3m(rad2deg([lat21_temp  ,lat21_temp ]), rad2deg([lon21_temp,  lon21_temp]),[0,alt21_temp/1000000],'color',[0.3 0.3 0.3]);
end
for i = 0:20/(time22(end)-time22(1)):1
    time22_temp = (time22(end)-time22(1))*i + time22(1);
    alt22_temp = interp1(time22,alt22,time22_temp);
    lat22_temp = interp1(time22,lat22,time22_temp);
    lon22_temp = interp1(time22,lon22,time22_temp);
    plot3m(rad2deg([lat22_temp  ,lat22_temp ]), rad2deg([lon22_temp,  lon22_temp]),[0,alt22_temp/1000000],'color',[0.3 0.3 0.3]);
end
for i = 0:20/(time3(end)-time3(1)):1
    time3_temp = (time3(end)-time3(1))*i + time3(1);
    alt3_temp = interp1(time3,alt3,time3_temp);
    lat3_temp = interp1(time3,lat3,time3_temp);
    lon3_temp = interp1(time3,lon3,time3_temp);
    plot3m(rad2deg([lat3_temp  ,lat3_temp ]), rad2deg([lon3_temp,  lon3_temp]),[0,alt3_temp/1000000],'color',[0.3 0.3 0.3]);
end
for i = 0:20/(timeexo(end)-timeexo(1)):1
    timeexo_temp = (timeexo(end)-timeexo(1))*i + timeexo(1);
    altexo_temp = interp1(timeexo,altexo,timeexo_temp);
    latexo_temp = interp1(timeexo,latexo,timeexo_temp);
    lonexo_temp = interp1(timeexo,lonexo,timeexo_temp);
    plot3m(rad2deg([latexo_temp  ,latexo_temp ]), rad2deg([lonexo_temp,  lonexo_temp]),[0,altexo_temp/1000000],'color',[0.3 0.3 0.3]);
end

zoom(30) 
%%

% MList = [5:0.1:10];
altlist = 20000:100:40000 ;
alphalist = 0:.1:6;

[alphagrid,altgrid] = ndgrid(alphalist,altlist);

v_temp = 2200;

for i = 1:numel(alphagrid)
        I = cell(1, ndims(alphagrid)); 
    [I{:}] = ind2sub(size(alphagrid),i);
    alpha_temp = alphagrid(I{1},I{2});
    alt_temp = altgrid(I{1},I{2});
    
    
[~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~] = VehicleModelCombined(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);

end

[~,~,~,~,~,~, ~, ~, Fd, ~,L,Fueldt,T,Isp,~,~,~,~] = VehicleModelCombined(0, 34000, 1600,auxdata,0,0,0,deg2rad(1),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0)
figure()
contourf(alphagrid,altgrid,(Tgrid-Fdgrid)./Fueldtgrid/9.81,1000,'LineWidth',0.)
title('Net Isp')
c = colorbar;

figure()
contourf(alphagrid,altgrid,Fdgrid,1000,'LineWidth',0.)
c = colorbar;

figure()
contourf(alphagrid,altgrid,Tgrid,1000,'LineWidth',0.)
title('Thrust')
c = colorbar;

figure()
contourf(alphagrid,altgrid,Fueldtgrid,1000,'LineWidth',0.)
c = colorbar;





%% SAVE FIGS

saveas(figure(311),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('ThirdStage',num2str(j),'.fig')]);
% print(figure(311),'ThirdStage','-dpng');
% movefile('ThirdStage.png',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(211),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('SecondStage',num2str(j),'.fig')]);
saveas(figure(221),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Return',num2str(j),'.fig')]);    
saveas(figure(2410),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Hamiltonian',num2str(j),'.fig')]);
saveas(figure(2420),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Validation',num2str(j),'.fig')]);
saveas(figure(212),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Forward1',num2str(j),'.fig')]);
saveas(figure(213),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Forward2',num2str(j),'.fig')]);
% saveas(figure(2100),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'eq.fig']);
% saveas(figure(2110),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ISP.fig']);
saveas(figure(2301),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('GroundTrack',num2str(j),'.fig')]);


close all



end




end



% Plot Cart Solution
addpath('.\axlabel')

% figure(401)
% surf(MList_EngineOff,AOAList_EngineOff,Cl_Grid_EngineOff,'facecolor','interp')
% xlabel('Mach no.')
% ylabel('Angle of Attack (deg)')
% zlabel('Cl')
% view(-20,30)
% align_axislabel([], gca)
% 
% 
% figure(402)
% surf(MList_EngineOff,AOAList_EngineOff,Cd_Grid_EngineOff,'facecolor','interp')
% view(-20,40)
% xlabel('Mach no.')
% ylabel('Angle of Attack (deg)')
% zlabel('Cd')
% view(-20,30)
% align_axislabel([], gca)

%------------------------------%
% Extract Solution from Output %
%------------------------------%
addpath('addaxis')

solution = output.result.solution;
time = solution.phase(1).time;
% altitude  = (solution.phase(1).state(:,1)-auxdata.Re);
% lon = solution.phase(1).state(:,2)*180/pi;
% lat  = solution.phase(1).state(:,3)*180/pi;
% speed     = solution.phase(1).state(:,4);
% fpa       = solution.phase(1).state(:,5)*180/pi;
% azimuth   = solution.phase(1).state(:,6)*180/pi;
% aoa       = solution.phase(1).state(:,7)*180/pi;
% bank      = solution.phase(1).state(:,8)*180/pi;
% mfuel      = solution.phase(1).state(:,9);
% throttle      = solution.phase(1).state(:,10);



c = ppval(auxdata.interp.c_spline,altitude); % Calculate speed of sound using atmospheric data
M = v./c; % Calculating Mach No (Descaled)

T0 = ppval(auxdata.interp.T0_spline, altitude); 

P0 = ppval(auxdata.interp.P0_spline, altitude);


[Isp,Fueldt,eq] = RESTint(M, aoa, auxdata,T0,P0);

Isp(M<5.0)=0;

rho = ppval(auxdata.interp.rho_spline,altitude); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

% Cd = auxdata.interp.Cd_spline(M,aoa);
% Cl = auxdata.interp.Cl_spline(M,aoa);
% %%%% Compute the drag and lift:
% A = auxdata.Stage2.Aref; % reference area (m^2)
% D = 0.5*Cd.*A.*rho.*speed.^2;
% L = 0.5*Cl.*A.*rho.*speed.^2;
%---------------%
% Plot Solution %
%---------------%

% figure('units','normalized','outerposition',[0.1 0.1 .9 .9])
% figure(1)
figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,1)
hold on


 plot(time,altitude/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
xlabel('time (s)');
addaxis(time,v/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Velocity (km/s)');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');

addaxis(time,azi,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Heading Angle (Deg)');


legend(  'Altitude', 'Velocity', 'Heading Angle', 'location', 'best');

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,2)
hold on
plot(time,rad2deg(aoa),'-','color','k', 'linewidth', 1.);
ylabel('Angle of Attack (deg)');
xlabel('time (s)')
throttle(M<5.0)=0;
addaxis(time,throttle*100,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Throttle (%)');

% addaxis(time,mfuel,'-.','color','k', 'linewidth', 1.);
% addaxislabel(3,'Fuel Mass (kg)');

addaxis(time,rad2deg(bank),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Bank Angle (Deg)');

addaxis(time,flapdeflection,'-.','color','k', 'linewidth', 1.2);
addaxislabel(4,'Flap Deflection (Deg)');
legend(  'Angle of Attack', 'Throttle' , 'Bank Angle', 'Flap Deflection');


figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
hold on
% subplot(3,1,3)
plot(time,M,'-','color','k', 'linewidth', 1.);
ylabel('Mach no.')
xlabel('time (s)')
addaxis(time,Isp,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Specific Impulse (s)');

addaxis(time,q,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

% addaxis(time,L./D,':','color','k', 'linewidth', 1.);
% addaxislabel(4,'L/D');

legend(  'Mach no.', 'Isp (potential)', 'q' );

print -depsc2 Trajectory.eps
print -dpng Trajectory.png


for i=1:length(output.meshhistory);
  mesh(i).points = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).iteration = i*ones(size(mesh(i).points));
end;

figure(4)
hold on

geoshow('landareas.shp','FaceColor',[0.5 .8 0.5])
xlim([min(lon)-3,max(lon)+3]);
ylim([min(lat)-3,max(lat)+3]);

plot(lon,lat,'-', 'color','k', 'linewidth', 1.5);

plot([lon(end) lon(1)],[lat(end) lat(1)],':', 'color','k', 'linewidth', 1.5);

txt1 = '2nd \rightarrow 3rd Stage Separation';
text(lon(1)-0.5,lat(1)-0.3,txt1)

txt2 = '2nd Stage Landing';
text(lon(end)+0.1,lat(end),txt2)

txt3 = 'Ascent';
text((lon(1)+lon(end))/2+0.1,(lat(1)+lat(end))/2,txt3)

txt3 = 'Return';
text(143,-10,txt3)

print -depsc2 rlvEntryLonLat.eps
print -dpng rlvEntryLonLat.png



figure(6);
for i=1:length(mesh);
  pp = plot(mesh(i).points,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location (Fraction of Interval)');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16,'FontName','Times');
grid on;

print -depsc2 MeshHistory.eps
print -dpng MeshHistory.png
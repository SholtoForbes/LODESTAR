% multi plotter
addpath('.\axlabel')
addpath('addaxis')
addpath('..\')

Atmosphere = dlmread('atmosphere.txt');
auxdata.Atmosphere = Atmosphere; 
auxdata.interp.c_spline = spline( auxdata.Atmosphere(:,1),  auxdata.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data


a = load('Solution-const6.mat');
b = load('Solution-const7.mat');

mFuel1  = a.solution.phase.state(1,9)
mFuel2  = b.solution.phase.state(1,9)

Re   = 6371203.92;                     % Equatorial Radius of Earth (m)

rad1  = a.solution.phase.state(:,1);
lon1  = rad2deg(a.solution.phase.state(:,2));
lat1  = rad2deg(a.solution.phase.state(:,3));
v1    = a.solution.phase.state(:,4);
fpa1  = a.solution.phase.state(:,5);
azi1  = a.solution.phase.state(:,6);
throttle1      = a.solution.phase.state(:,10);
bank1      = a.solution.phase.state(:,8);
t1 = a.solution.phase.time;

alt1 = rad1 - Re;

c1 = ppval(auxdata.interp.c_spline,alt1); % Calculate speed of sound using atmospheric data
M1 = v1./c1; % Calculating Mach No (Descaled)
throttle1(M1<5.1)=0;

rad2  = b.solution.phase.state(:,1);
lon2  = rad2deg(b.solution.phase.state(:,2));
lat2  = rad2deg(b.solution.phase.state(:,3));
v2    = b.solution.phase.state(:,4);
fpa2  = b.solution.phase.state(:,5);
azi2  = b.solution.phase.state(:,6);
throttle2      = b.solution.phase.state(:,10);
bank2      = b.solution.phase.state(:,8);
t2 = b.solution.phase.time;

alt2 = rad2 - Re;

c2 = ppval(auxdata.interp.c_spline,alt2); % Calculate speed of sound using atmospheric data
M2 = v2./c2; % Calculating Mach No (Descaled)
throttle2(M2<5.1)=0;

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
hold on
plot(t1,(rad1-Re)/1000,'-','color','k', 'linewidth', 1.);
plot(t2,(rad2-Re)/1000,'--','color','k', 'linewidth', 1.);
xlabel('Time (s)');
ylabel('Altitude (km)');

addaxis(t1,v1,'-','color',[0.5 0 0], 'linewidth', 1.);
addaxisplot(t2,v2,2,'--','color',[0.5 0 0], 'linewidth', 1.);
addaxislabel(2,'Velocity(m/s)');

addaxis(t1,throttle1,'-','color',[0 0 0.5], 'linewidth', 1.);
addaxisplot(t2,throttle2,3,'--','color',[0 0 0.5], 'linewidth', 1.);
addaxislabel(3,'Throttle');

% legend( '110% Cd Alt','90% Cd Alt','110% Cd v','90% Cd v','110% Cd Throttle','90% Cd Throttle')
legend( '110% Isp Alt','90% Isp Alt','110% Isp v','90% Isp v','110% Isp Throttle','90% Isp Throttle')
set(gca,'box','off')
figure(3)
hold on
geoshow('landareas.shp','FaceColor',[0.5 .8 0.5])
xlim([min(lon1)-3,max(lon1)+3]);
ylim([min(lat1)-3,max(lat1)+3]);

plot(lon1,lat1,'-', 'color','k', 'linewidth', 1.5);
plot(lon2,lat2,'--', 'color','k', 'linewidth', 1.5);

txt1 = '2nd \rightarrow 3rd Stage Separation';
text(lon1(1)-0.5,lat1(1)-0.3,txt1)

txt2 = '2nd Stage Landing';
text(lon1(end)+0.1,lat1(end),txt2)


figure()
hold on
plot(t1,bank1)
plot(t2,bank2)

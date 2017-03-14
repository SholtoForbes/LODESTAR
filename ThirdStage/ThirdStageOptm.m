function mpayload = ThirdStageOptm(k,j,u, phi0, zeta0)

mScale = 1; % This needs to be manually changed in altitude and velocity files as well
% x0 = [1200*mScale deg2rad(15) deg2rad(15)] % 
% x0 = [1500  deg2rad(13) deg2rad(13)];

[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 0 10],k,j,u, phi0, zeta0);
% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 0 0 10],k,j,u, phi0, zeta0);
AoA_max

% x0 = [1600  deg2rad(10) deg2rad(10) deg2rad(10) deg2rad(10) deg2rad(10) 100];
% x0 = [1600  deg2rad(10) deg2rad(10) deg2rad(10) deg2rad(10) 200];
x0 = [2550/10000  AoA_max AoA_max AoA_max AoA_max 250/1000];

% x0 = [2700  deg2rad(7.5) deg2rad(10)]
% x0 = [2750  deg2rad(6) AoA_max AoA_max deg2rad(11) deg2rad(10)]
% x0 = [2700  deg2rad(8) deg2rad(10) deg2rad(10) deg2rad(10)]

% x0 = [2750  deg2rad(14) deg2rad(10) deg2rad(14) deg2rad(12)]
% x0 = [2750  deg2rad(6) AoA_max-0.01 AoA_max-0.01 deg2rad(8)]


% x0 = [2750  (deg2rad(14)+(deg2rad(6)-deg2rad(14))*j/0.05) (deg2rad(10)+((AoA_max-0.01)-deg2rad(10))*j/0.05) (deg2rad(14)+((AoA_max-0.01)-deg2rad(14))*j/0.05) (deg2rad(12)+(deg2rad(8)-deg2rad(12))*j/0.05)]
% x0 = [2750  (deg2rad(14)+(deg2rad(6)-deg2rad(14))*j/0.05) (deg2rad(10)+((AoA_max-0.01)-deg2rad(10))*j/0.05) (deg2rad(14)+((AoA_max-0.01)-deg2rad(14))*j/0.05) (deg2rad(12)+(deg2rad(8)-deg2rad(12))*j/0.05) 200]


% x0 = [1600  AoA_max];
% x0 = [1650 (AoA_max-0.01)*10000];
% x0 = [1650];
options.Display = 'iter-detailed';
% options.Algorithm = 'sqp';
options.TypicalX = x0;
% options.UseParallel = 1;
% options.Algorithm = 'active-set';

% options.TolFun = 1e-10;
% options.TolX = 1e-10;

% k = 35500;
% j = 0.05;
% u = 2840;
% x = fminsearch(@(x)Payload(x,k,j,u, phi0, zeta0),x0,options);
% x = fminunc(@(x)Payload(x,k,j,u, phi0, zeta0),x0,options);
%  x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[1400 0 0 0 0 0 50],[2900 AoA_max AoA_max AoA_max AoA_max AoA_max 300],@(x)Constraint(x,k,j,u, phi0, zeta0),options);
% x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[1400 0 0 0 0 50],[3200 AoA_max AoA_max AoA_max AoA_max 250],@(x)Constraint(x,k,j,u, phi0, zeta0),options);
x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[2200/10000 0 0 0 0 50/1000],[3000/10000 AoA_max AoA_max AoA_max AoA_max 300/1000],@(x)Constraint(x,k,j,u, phi0, zeta0),options);

mfuel_burn = x(1)
AoA_control1 = x(2)
% AoA_control2 = x(3)


[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x,k,j,u, phi0, zeta0);
mpayload
zeta(end)
figure(9)
xlabel('time (s)')
set(gcf,'position',[300 300 800 600])


subplot(2,1,1);
hold on
plot(t, Alt/100, 'LineStyle', '-','Color','k', 'lineWidth', 1.3)
plot(t,v, 'LineStyle', '--','Color','k', 'lineWidth', 1.2)
plot(t, m, 'LineStyle', ':','Color','k', 'lineWidth', 1.4)
legend(  'Altitude (km x 100)',  'Mass (kg)', 'Velocity (m/s)');
subplot(2,1,2);
hold on
plot(t, rad2deg(gamma), 'LineStyle', '--','Color','k', 'lineWidth', 1.3)
plot(t(1:end-1),q/10000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
plot(t(1:end-1),rad2deg(Alpha)/10, 'LineStyle', '-','Color','k', 'lineWidth', 1.1)
legend(  'Trajectory Angle (degrees)','Dynamic Pressure (kPa) x 10','Angle of Attack (deg) x 10');

% legend(  'Altitude (km x 100)', 'Trajectory Angle (degrees)', 'Mass (kg x 10^3)', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (kPa) x 10','Angle of Attack (deg) x 10');
ylim([0 8])
xlim([0 t(end)])

dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')

Integrated_Drag = cumtrapz(t(1:end-1),D) ;
Integrated_Drag(end)
end
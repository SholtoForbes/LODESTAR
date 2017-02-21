function mpayload = ThirdStageOptm(k,j,u, phi0, zeta0)

mScale = 1; % This needs to be manually changed in altitude and velocity files as well
% x0 = [1200*mScale deg2rad(15) deg2rad(15)] % 
% x0 = [1500  deg2rad(13) deg2rad(13)];

[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0],k,j,u, phi0, zeta0);

AoA_max

% x0 = [1600  AoA_max AoA_max];
x0 = [1600  AoA_max];
% x0 = [1650 (AoA_max-0.01)*10000];
% x0 = [1650];
options.Display = 'iter';
% options.Algorithm = 'sqp';
options.TypicalX = x0;
% options.Algorithm = 'active-set';

options.TolFun = 1;
options.TolX = 1;

% k = 35500;
% j = 0.05;
% u = 2840;
x = fminsearch(@(x)Payload(x,k,j,u, phi0, zeta0),x0,options);
% x = fminunc(@(x)Payload(x,k,j,u, phi0, zeta0),x0,options);
% x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[1100 0],[1800 AoA_max*10000],[],options);
% x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[1100],[1800],[],options);
mfuel_burn = x(1)
AoA_control1 = x(2)/10000
% AoA_control2 = x(3)


[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x,k,j,u, phi0, zeta0);
mpayload
zeta(end)
figure(9)
xlabel('time (s)')
set(gcf,'position',[300 300 800 600])
hold on
plot(t, Alt/10000, 'LineStyle', '-','Color','k', 'lineWidth', 1.3)
plot(t, rad2deg(gamma), 'LineStyle', '--','Color','k', 'lineWidth', 1.3)
plot(t, m/100, 'LineStyle', ':','Color','k', 'lineWidth', 1.4)
plot(t,v/100, 'LineStyle', '--','Color','k', 'lineWidth', 1.2)
plot(t(1:end-1),q/1000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
plot(t(1:end-1),rad2deg(Alpha), 'LineStyle', '-','Color','k', 'lineWidth', 1.1)

legend(  'Altitude (km x 10)', 'Trajectory Angle (degrees)', 'Mass (kg x 10^2)', 'Velocity (m/s x 10^2)', 'Dynamic Pressure (kPa)','Angle of Attack (deg)');
ylim([0 60])
xlim([0 t(end)])

dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')

Integrated_Drag = cumtrapz(t(1:end-1),D) ;
Integrated_Drag(end)
end
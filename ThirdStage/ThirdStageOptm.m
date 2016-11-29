function mpayload = ThirdStageOptm(k,j,u, phi0, zeta0)

mScale = 1; % This needs to be manually changed in altitude and velocity files as well
% x0 = [1200*mScale deg2rad(15) deg2rad(15)] % 
% x0 = [1500  deg2rad(13) deg2rad(13)];
x0 = [1500  deg2rad(13)];
options.Display = 'iter';
options.TolFun = 0.1;
options.TolX = 0.1;

% k = 35500;
% j = 0.05;
% u = 2840;
x = fminsearch(@(x)Payload(x,k,j,u, phi0, zeta0),x0,options);
% x = fmincon(@(x)Payload(x,k,j,u),x0,[1],[1800],[],[],[],[],[],options);

mfuel_burn = x(1)
AoA_control = [x(2)]

[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma] = ThirdStageSim(x,k,j,u, phi0, zeta0);


figure(5)
xlabel('time (s)')
set(gcf,'position',[300 300 800 600])
hold on
plot(t, Alt/1000, 'LineStyle', '-','Color','k', 'lineWidth', 2.0)
plot(t, rad2deg(gamma), 'LineStyle', '--','Color','k', 'lineWidth', 1.0)
plot(t, m/100, 'LineStyle', ':','Color','k', 'lineWidth', 2.0)
plot(t,v/100, 'LineStyle', '--','Color','k', 'lineWidth', 2.0)
plot(t(1:end-1),q/1000, 'LineStyle', '-.','Color','k', 'lineWidth', 2.0)
% plot(t(1:end-1),rad2deg(Alpha), 'LineStyle', '-','Color','k', 'lineWidth', 1.0)

legend(  'Altitude (km)', 'Trajectory Angle (degrees)', 'Mass (kg x 10^2)', 'Velocity (m/s x 10^2)', 'Dynamic Pressure (kPa)')

end
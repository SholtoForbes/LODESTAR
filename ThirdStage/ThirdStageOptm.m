function [mpayload, x] = ThirdStageOptm(k,j,u, phi0, zeta0)

mScale = 1; % This needs to be manually changed in altitude and velocity files as well


% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 20],k,j,u, phi0, zeta0);


[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([10 0 0 0 0 20],k,j,u, phi0, zeta0);
AoA_max


% x0 = [2590/10000  AoA_max*ones(1,16) 250/1000]; % Works for 36km 

% x0 = [2590/10000  AoA_max*ones(1,16)-deg2rad(.5) 250/1000]; % this works well for 33 and 34km, except for 0 gamma 34km, doesnt work for 36km
% works for 35km 0 gamma

% 
% x0 = [2590/10000  AoA_max*ones(1,16)-deg2rad(2) 250/1000]; % seems to do well for 0 gamma 33km and 34km

% x0 = [2590/10000  AoA_max*ones(1,16)-deg2rad(3.5) 250/1000]; % 

% this works pretty well, just not at 33km 0 gamma yet
% nodesalt = [33000; 33000; 34000; 36000 ;36000];
% nodesgam = [0;0.05; 0; 0; 0.05];
% vals =  [deg2rad(2);deg2rad(.5); deg2rad(2); 0 ;0];
% interp = scatteredInterpolant(nodesalt,nodesgam,vals);
% x0 = [2590/10000  AoA_max*ones(1,16)-interp(k,j) 250/1000]; % this problem is extremely sensitive to initital guess! mostly at low altitude low gamma


%burn fuel and end aoa time are currently contrained in thirdstagesim


% x0 = [2590/10000  AoA_max*ones(1,5)  0*ones(1,1) AoA_max*ones(1,10)  250/1000];

options.Display = 'iter-detailed';
options.Algorithm = 'sqp';
% options.DiffMinChange = 0.0005;
% options.TypicalX = x0;
% options.UseParallel = 1;
% options.Algorithm = 'active-set';

options.TolFun = 1e-3;
options.TolX = 1e-3;

mpayload = 0;
x=0;

for i = 2
% x0 = [2590/10000  AoA_max*ones(1,16)-deg2rad(i/2) 250/1000];
% x0 = [2590/10000  AoA_max*ones(1,16) 250/1000];
x0 = [2590/10000  AoA_max*ones(1,30) 250/1000];

% x0 = [2590/10000  AoA_max*ones(1,10) 300000/1e7 250/1000];

% x0 = [2590/10000  0*ones(1,30) AoA_max-deg2rad(1)];
% options.DiffMinChange = 0.0006 + 0.0001*i;
options.DiffMinChange =  0.001*i;
[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[2300/10000 deg2rad(0)*ones(1,30) 200/1000],[3000/10000 AoA_max*ones(1,30) 300/1000],@(x)Constraint(x,k,j,u, phi0, zeta0),options);

% [x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[2300/10000 deg2rad(0)*ones(1,10) 160000/1e7 200/1000],[3000/10000 AoA_max*ones(1,10) 580000/1e7 300/1000],@(x)Constraint(x,k,j,u, phi0, zeta0),options);

% [x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[2300/10000 deg2rad(-1)*ones(1,30) 0],[3000/10000 deg2rad(1)*ones(1,30) AoA_max],@(x)Constraint(x,k,j,u, phi0, zeta0),options);
exitflag
[AltF, vF, Alt, v, t, mpayload_temp, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x_temp,k,j,u, phi0, zeta0);

if mpayload_temp > mpayload && (exitflag ==1 || exitflag ==2|| exitflag ==3)
    mpayload = mpayload_temp;
    x = x_temp;
end
end
%     mpayload = mpayload_temp;
%     x = x_temp;
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x,k,j,u, phi0, zeta0);

% x = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[2200/10000 deg2rad(0)*ones(1,16) 200/1000],[3000/10000 AoA_max*ones(1,16) 270/1000],@(x)Constraint(x,k,j,u, phi0, zeta0),options);
% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x,k,j,u, phi0, zeta0);

mfuel_burn = x(1)
AoA_control1 = x(2)
% AoA_control2 = x(3)
x(end)



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
% plot(t(1:end),rad2deg(Alpha)/10, 'LineStyle', '-','Color','k', 'lineWidth', 1.1)
legend(  'Trajectory Angle (degrees)','Dynamic Pressure (kPa) x 10','Angle of Attack (deg) x 10');

% legend(  'Altitude (km x 100)', 'Trajectory Angle (degrees)', 'Mass (kg x 10^3)', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (kPa) x 10','Angle of Attack (deg) x 10');
ylim([0 8])
xlim([0 t(end)])

dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')

Integrated_Drag = cumtrapz(t(1:end-1),D) ;
Integrated_Drag(end)
end
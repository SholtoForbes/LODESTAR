function [mpayload, x, zeta, phi,Alt,v,t,Alpha,m,gamma,q] = ThirdStageOptm(k,j,u, phi0, zeta0)
% A function to calculate the optimal trajectory of the SPARTAN third stage.
% Sholto Forbes-Spyratos

% The settings in this function should match those in the MatrixParallel routine

mScale = 1; % This needs to be manually changed in altitude and velocity files as well

% This is used to calculate the max AoA, not the best coding here
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 20],k,j,u, phi0, zeta0);
AoA_max


% options.Display = 'iter-detailed';
options.Algorithm = 'sqp';
options.FunValCheck = 'on';
% options.ScaleProblem = 'obj-and-constr'
% options.DiffMinChange = 0.0005;
% options.TypicalX = x0;
% options.UseParallel = 1;
% options.Algorithm = 'active-set';


options.TolFun = 1e-3;
options.TolX = 1e-3;

mpayload = 0;
x=0;
% Run each case for a range of initial guesses and DiffMinChange values.
% This mostly ensures that the optimal solution will be found.
for i3 = 0:.5:6
for i2 = 0:10

i4=0;
x0 = [AoA_max*ones(1,10)-i4*AoA_max*0.01 250/10000+i2*5/10000]; 

% Initiate optimiser
options.DiffMinChange = 0.0005*i3;
[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0),x0,[],[],[],[],[deg2rad(0)*ones(1,10) 200/10000],[AoA_max*ones(1,10) 350/10000],@(x)Constraint(x,k,j,u, phi0, zeta0),options);

exitflag
[AltF, vF, Alt, v, t, mpayload_temp, Alpha, m,AoA,q,gamma,D,AoA_max,zeta,inc] = ThirdStageSim(x_temp,k,j,u, phi0, zeta0);

if mpayload_temp > mpayload && (exitflag ==1 || exitflag ==2|| exitflag ==3)
    mpayload = mpayload_temp;
    x = x_temp;
end
end
end

[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max,zeta,phi,inc] = ThirdStageSim(x,k,j,u, phi0, zeta0);


mfuel_burn = x(1)
AoA_control1 = x(2)
x(end)



mpayload
zeta(end)
figure(301)
xlabel('time (s)')
set(gcf,'position',[300 300 800 600])

%% Plotting
subplot(2,1,1);
hold on
plot(t, Alt/100, 'LineStyle', '-','Color','k', 'lineWidth', 1.3)
plot(t,v, 'LineStyle', '--','Color','k', 'lineWidth', 1.2)
plot(t, m, 'LineStyle', ':','Color','k', 'lineWidth', 1.4)
legend(  'Altitude (km x 10)', 'Velocity (m/s)',  'Mass (kg)');
subplot(2,1,2);
hold on
plot(t, rad2deg(gamma), 'LineStyle', '--','Color','k', 'lineWidth', 1.3)
plot(t(1:end-1),q/10000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
plot(t(1:end-1),rad2deg(Alpha)/10, 'LineStyle', '-','Color','k', 'lineWidth', 1.1)
legend(  'Trajectory Angle (degrees)','Dynamic Pressure (kPa) x 10','Angle of Attack (deg) x 10');
ylabel('Time (s)');
ylim([-1 8])
xlim([0 t(end)])
% Write data to file
dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')

Integrated_Drag = cumtrapz(t(1:end-1),D) ;
Integrated_Drag(end)
end
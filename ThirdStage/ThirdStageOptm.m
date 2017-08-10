function [mpayload, x, zeta, phi,Alt,v,t,Alpha,m,gamma,q,Vec_angle,T,CL,L,AltF_actual,vF] = ThirdStageOptm(k,j,u, phi0, zeta0)
% A function to calculate the optimal trajectory of the SPARTAN third stage.
% Sholto Forbes-Spyratos

% The settings in this function should match those in the MatrixParallel routine

mScale = 1; % This needs to be manually changed in altitude and velocity files as well

% This is used to calculate the max AoA, not the best coding here
% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 20],k,j,u, phi0, zeta0);
% AoA_max(1)


% options.Display = 'iter-detailed';
options.Algorithm = 'sqp';
% options.FunValCheck = 'on';
options.ScaleProblem = 'obj-and-constr'
% options.DiffMinChange = 0.0005;
% options.TypicalX = x0;
% options.UseParallel = 1;
% options.Algorithm = 'active-set';


% options.TolFun = 1e-4;
% options.TolX = 1e-4;

mpayload = 0;
x=0;
% Run each case for a range of initial guesses and DiffMinChange values.
% This mostly ensures that the optimal solution will be found.

% count = 1

    
    for i5 = 0;
    for i2 = 0:.25:1.5;
        for i4 = 0;

i3 = 1;
% for i4 = 0:5;
%     for i3 = 0:2;
%     for i5 = 0;

[i5 i2 i4 ]
        
% count = count+1
AoA_max_abs = deg2rad(15); % maximum angle of attack


num_div = 20+i5; 
% num_div = 15-i4;

% x0 = [deg2rad(11)*ones(1,num_div)+deg2rad(i4) 2650/10000 245/1000];
x0 = [deg2rad(14)*ones(1,num_div)+deg2rad(i4) 2450/10000+i2*50/10000 220/1000];

% x0 = [deg2rad(11)*ones(1,num_div)+deg2rad(i4) 2650/10000 85/1000];


% x0 = [(30+i4*10)/1000 deg2rad(20)*ones(1,num_div) 2800/10000 230/1000];
% x0 = [deg2rad(10)*ones(1,num_div)+deg2rad(i4) 2700/10000];

% Initiate optimiser
options.DiffMinChange = 0.0005*i3;

lb = [deg2rad(0)*ones(1,num_div) 2400/10000  200/1000];
ub = [AoA_max_abs*ones(1,num_div) 2900/10000  240/1000];
% lb = [deg2rad(0)*ones(1,num_div) 2550/10000  70/1000];
% ub = [AoA_max_abs*ones(1,num_div) 2900/10000  90/1000];


% lb = [20/1000 deg2rad(0)*ones(1,num_div) 2500/10000  200/1000];
% ub = [150/1000 AoA_max_abs*ones(1,num_div) 2900/10000  250/1000];
% lb = [deg2rad(0)*ones(1,num_div) 2400/10000];
% ub = [AoA_max_abs*ones(1,num_div) 2900/10000 ];


[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0,lb,num_div),x0,[],[],[],[],lb,ub,@(x)Constraint(x,k,j,u, phi0, zeta0,lb,num_div),options);

 
exitflag
[AltF_actual, vF, Alt, v, t, mpayload_temp, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L] = ThirdStageSim(x_temp,k,j,u, phi0, zeta0, lb,num_div);
mpayload_temp
AltF_actual
% Alpha
Vec_angle_constraint = max(Vec_angle - deg2rad(7)); % check for thrust vector validity (constraints are not necessarily satisfied)

if mpayload_temp > mpayload && (exitflag ==1 || exitflag ==2|| exitflag ==3) && Vec_angle_constraint <= 0 
    mpayload = mpayload_temp;
    x = x_temp;
    num_div_best = num_div;
end
x
mpayload
end
end
    end
%     end
x
[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div_best);


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
plot(t(1:end-1),rad2deg(Vec_angle)/10, 'LineStyle', ':','Color','k', 'lineWidth', 1.1)
legend(  'Trajectory Angle (degrees)','Dynamic Pressure (kPa) x 10','Angle of Attack (deg)x10', 'Thrust Vector Angle (deg)x10');
ylabel('Time (s)');
ylim([-1 8])
xlim([0 t(end)])
% Write data to file
dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')

Integrated_Drag = cumtrapz(t(1:end-1),D) ;
Integrated_Drag(end)
end
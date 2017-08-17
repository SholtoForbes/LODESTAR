function [mpayload, x, zeta, phi,Alt,v,t,Alpha,m,gamma,q,Vec_angle,T,CL,L,AltF_actual,vF,inc_diff] = ThirdStageOptm(k,j,u, phi0, zeta0)
% A function to calculate the optimal trajectory of the SPARTAN third stage.
% Sholto Forbes-Spyratos

% The settings in this function should match those in the MatrixParallel routine

mScale = 1; % This needs to be manually changed in altitude and velocity files as well

% This is used to calculate the max AoA, not the best coding here
% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 20],k,j,u, phi0, zeta0);
% AoA_max(1)
temp_payload = [];
mpayload_mat = [];
x_mat = {};
num_div_mat = [];
options.Display = 'none';
options.Algorithm = 'sqp';
% options.FunValCheck = 'on';
options.ScaleProblem = 'obj-and-constr'
% options.DiffMinChange = 0.0005;
% options.TypicalX = x0;
% options.UseParallel = 1;
% options.Algorithm = 'active-set';


options.TolFun = 1e-4;
options.TolX = 1e-4;
options.DiffMinChange = 0.0005;
%  options.DiffMinChange = 0.01;
mpayload = 0;
x=0;
% Run each case for a range of initial guesses and DiffMinChange values.
% This mostly ensures that the optimal solution will be found.

successflag = 0; % flags if any runs are successful


% count = 1
AoA_max_abs = deg2rad(20); % maximum angle of attack




    parfor i = 1:18
%         i4=0;
        i5=0;
        
        if i <= 6
i3 = i;
i4 = 0;

        elseif i <=12
 i3 = i-7;
i4 = 1;  
        elseif i <=18
             i3 = i-13;
i4 = 2;  
        elseif i <= 24
                      i3 = i-19;
i4 = 3;     
        elseif i <= 32
                      i3 = i-25;
i4 = 4;     
        end

num_div = 20+i5;

%     lb = [deg2rad(0)*ones(1,num_div) 3100/10000 2400/10000  200/1000];
% ub = [AoA_max_abs*ones(1,num_div) 3300/10000 2700/10000  240/1000];
% 
% x0 = [deg2rad(11)*ones(1,num_div)+deg2rad(i4) 3150/10000 2425/10000+i3*25/10000 220/1000];

%     lb = [deg2rad(0)*ones(1,num_div) 2400/10000  200/1000];
% ub = [AoA_max_abs*ones(1,num_div) 2900/10000  240/1000];
% 
% x0 = [deg2rad(12)*ones(1,num_div)+deg2rad(i4) 2425/10000+i3*25/10000 220/1000];


 lb = [deg2rad(0)*ones(1,num_div) 1300/10000  70/1000];
ub = [AoA_max_abs*ones(1,num_div) 2200/10000  130/1000];

x0 = [deg2rad(18)*ones(1,num_div)+deg2rad(i4) 1900/10000+i3*25/10000 110/1000];


% Initiate optimiser


[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u, phi0, zeta0,lb,num_div),x0,[],[],[],[],lb,ub,@(x)Constraint(x,k,j,u, phi0, zeta0,lb,num_div),options);

 [~, ~,~, ~, ~, mpayload_temp, ~, ~,~,~,~,~,~,~,~, ~,Vec_angle] = ThirdStageSim(x_temp,k,j,u, phi0, zeta0, lb,num_div);

% exitflag
% [AltF_actual, vF, Alt, v, t, mpayload_temp, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L] = ThirdStageSim(x_temp,k,j,u, phi0, zeta0, lb,num_div);

% mpayload_temp
% AltF_actual
% Alpha
Vec_angle_constraint = max(Vec_angle - deg2rad(8)); % check for thrust vector validity (constraints are not necessarily satisfied)

if mpayload_temp > mpayload && (exitflag ==1 || exitflag ==2|| exitflag ==3) && Vec_angle_constraint <= 0 
mpayload_mat(i) = mpayload_temp;
x_mat{i} = x_temp;
num_div_mat(i) = num_div;
successflag(i) = 1;
else
mpayload_mat(i) = 0;
x_mat{i} = 0*ones(1,length(x_temp));
num_div_mat(i) = num_div;
end

    end
%     end
[~,best_location] = max(mpayload_mat);
x = x_mat{best_location};
num_div_best = num_div_mat(best_location);
    lb = [deg2rad(0)*ones(1,num_div_best) 2400/10000  200/1000];
ub = [AoA_max_abs*ones(1,num_div_best) 2900/10000  240/1000];
x
% [AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div_best);
if any(successflag)
[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L,inc_diff] = ThirdStageSim(x,k,j,u, phi0, zeta0, lb,num_div_best);
else
AltF_actual=0;
    vF=0; 
    Alt=0; 
    v=0; 
    t=0; 
    mpayload=0; 
    Alpha=0;
    m=0;
    AoA_init=0;
    q=0;
    gamma=0;
    D=0;
AoA_max=0;
    zeta=0;
    phi=0; 
    
inc=0;
    Vec_angle=0;
    T=0;
    CL=0;
    L=0;
    inc_diff=0;
end

mfuel_burn = x(1)
AoA_control1 = x(2)
x(end)



mpayload
zeta(end)
figure(301)
xlabel('time (s)')
set(gcf,'position',[300 300 800 600])

%% Plotting
% subplot(2,1,1);
% hold on
% plot(t, Alt/100, 'LineStyle', '-','Color','k', 'lineWidth', 1.3)
% plot(t,v, 'LineStyle', '--','Color','k', 'lineWidth', 1.2)
% plot(t, m, 'LineStyle', ':','Color','k', 'lineWidth', 1.4)
% legend(  'Altitude (km x 10)', 'Velocity (m/s)',  'Mass (kg)');
% subplot(2,1,2);
% hold on
% 
% plot(t, rad2deg(gamma), 'LineStyle', '--','Color','k', 'lineWidth', 1.3)
% plot(t(1:end-1),q/10000, 'LineStyle', '-.','Color','k', 'lineWidth', 1.0)
% 
% plot(t(1:end-1),rad2deg(Alpha)/10, 'LineStyle', '-','Color','k', 'lineWidth', 1.1)
% plot(t(1:end-1),rad2deg(Vec_angle)/10, 'LineStyle', ':','Color','k', 'lineWidth', 1.1)
% legend(  'Trajectory Angle (degrees)','Dynamic Pressure (kPa) x 10','Angle of Attack (deg)x10', 'Thrust Vector Angle (deg)x10');
% ylabel('Time (s)');
% ylim([-1 8])
% xlim([0 t(end)])
% % Write data to file
% dlmwrite('ThirdStageData',[t.', Alt.', v.', m.',[q q(end)].',gamma.',[D D(end)].',zeta.'], ' ')
% 
% Integrated_Drag = cumtrapz(t(1:end-1),D) ;
% Integrated_Drag(end)
end
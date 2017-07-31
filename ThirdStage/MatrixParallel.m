%% Routine to calculate a matrix of payload mass to orbit
% Sholto Forbes-Spyratos

clear all
mat = [];

u = [2825:25:2925] % velocity range. The routine is parallelised around this velocity list
% u = [2825:25:3000] % velocity range. The routine is parallelised around this velocity list

phi0 = -0.13 % initial latitude, this has very minimal effect
% zeta0 = 1.69 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
 zeta0 = 1.76 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
for k = [33000:1000:38000] % altitude range

for j = [0 0.0125 0.025 0.0375 0.05] % trajectory angle range
% for j = [0.025 0.0375 0.05]
    
temp_guess_no = 1;

phi0
zeta0
k
j
AltF = [];
vF = [];
Alt = [];
t = [];
mpayload = [];
Alpha = [];
m = [];
AoA = [];
options = cell(1,8);

parfor i = 1:length(u)

% This is used to calculate the max AoA, not the best coding here
% [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0*ones(1,5) 10],k,j,u(i), phi0, zeta0);

mpayload(i) = 0;
options{i}.Display = 'none';
% options{i}.Display = 'iter';
options{i}.Algorithm = 'sqp';
% options{i}.Algorithm = 'active-set';
options{i}.ScaleProblem = 'obj-and-constr';

options{i}.TolFun = 1e-4;
options{i}.TolX = 1e-4;
    
% Run each case for a range of initial guesses and DiffMinChange values.
% This mostly ensures that the optimal solution will be found.
% for i3 = 0:.5:12
% for i2 = 0:2.5:10
%  for i2 = 0:1:10
    
% for i3 = 0:.25:6
% %     for i3 = 0:.1:8
% for i4 = 0:2

% for i3 = 0:.5:7
% for i2 = 0:2.5:5
% for i4 = 0:2;

% for i3 = 0:3
% for i4 = 0:.25:3;
    
    for i4 = 0:.25:1;
    for i3 = 0:2
%         i3 = 3;
   for i5 = 0:1;
% i2 = 1;
% i4 = 0;

AoA_max_abs = deg2rad(15);
% 
% x0 = [AoA_max(1)*ones(1,10)-i4*AoA_max(1)*0.01 250/10000+i2*5/10000]; % initial guess uses first max aoa
% 
% if AoA_max(1) > AoA_max_abs
%     x0 = [AoA_max_abs*ones(1,10)-i4*AoA_max(1)*0.01 250/10000+i2*5/10000];
% end

% x0 = [deg2rad(5)*ones(1,10)+deg2rad(i4) 250/10000]; % initial guess uses first max aoa
% x0 = [deg2rad(5)*ones(1,10) 230/10000];

% x0 = [AoA_max*ones(1,10)-i4*AoA_max*0.01];

num_div = 20+i5;
x0 = [deg2rad(14)*ones(1,num_div)+deg2rad(i4) 2800/10000 230/1000];


options{i}.DiffMinChange = 0.0005*i3;

lb = [deg2rad(0)*ones(1,num_div) 2500/10000  200/1000];

% Initiate optimiser
% [x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[deg2rad(0)*ones(1,10) 200/10000],[AoA_max_abs*ones(1,10) 350/10000],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options{i});
%  [x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[deg2rad(0)*ones(1,10)],[AoA_max_abs*ones(1,10)],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options{i});
% [x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[deg2rad(0)*ones(1,8)  200/10000],[AoA_max_abs*ones(1,8)  240/10000],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options{i});
[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0,lb,num_div),x0,[],[],[],[],lb,[AoA_max_abs*ones(1,num_div) 2900/10000  240/1000],@(x)Constraint(x,k,j,u(i), phi0, zeta0,lb,num_div),options{i});

[AltF(i), vF(i), Alt, v, t, mpayload_temp, Alpha, m,AoA,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle] = ThirdStageSim(x_temp,k,j,u(i), phi0, zeta0,lb,num_div);
Vec_angle_constraint = max(Vec_angle - deg2rad(20));

if mpayload_temp > mpayload(i) && (exitflag ==1 || exitflag ==2|| exitflag ==3) && Vec_angle_constraint <= 0 
    mpayload(i) = mpayload_temp; % if payload improves, set new max payload
end
end
    end
    end
temp_payload(i) = mpayload(i);
end
mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',temp_payload.']];
temp_guess_no = temp_guess_no + 1;
end
end

dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')
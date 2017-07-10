%% Routine to calculate a matrix of payload mass to orbit
% Sholto Forbes-Spyratos

clear all
mat = [];

u = [2850:25:2925] % velocity range. The routine is parallelised around this velocity list
% u = [2825:25:3000] % velocity range. The routine is parallelised around this velocity list

phi0 = -0.13 % initial latitude, this has very minimal effect
% zeta0 = 1.69 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
 zeta0 = 1.76 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
for k = [33000:1000:37000] % altitude range

for j = [0 0.0125 0.025 0.0375 0.05] % trajectory angle range

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
[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0*ones(1,5) 10],k,j,u(i), phi0, zeta0);

mpayload(i) = 0;
% options{i}.Display = 'final';
% options{i}.Display = 'iter';
options{i}.Algorithm = 'sqp';
% options{i}.Algorithm = 'active-set';
% options(i).ScaleProblem = 'iter-and-constr'

options{i}.TolFun = 1e-3;
options{i}.TolX = 1e-3;
    
% Run each case for a range of initial guesses and DiffMinChange values.
% This mostly ensures that the optimal solution will be found.
for i3 = 0:.5:6
for i2 = 0:10
% for i3 = 0:2:6
% for i2 = 0:2.5:10
    

i4 = 0;
x0 = [AoA_max*ones(1,10)-i4*AoA_max*0.01 250/10000+i2*5/10000]; % initial guess

options{i}.DiffMinChange = 0.0005*i3;

% Initiate optimiser
[x_temp,fval,exitflag] = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[deg2rad(0)*ones(1,10) 200/10000],[AoA_max*ones(1,10) 350/10000],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options{i});
[AltF(i), vF(i), Alt, v, t, mpayload_temp, Alpha, m,AoA,q,gamma,D,AoA_max,zeta] = ThirdStageSim(x_temp,k,j,u(i), phi0, zeta0);

if mpayload_temp > mpayload(i) && (exitflag ==1 || exitflag ==2|| exitflag ==3)
    mpayload(i) = mpayload_temp; % if payload improves, set new max payload
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
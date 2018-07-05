%% Routine to calculate a matrix of payload mass to orbit
% Sholto Forbes-Spyratos

clear all
mat = [];

% u = [2850:25:2925] % velocity range. The routine is parallelised around this velocity list
u = [2850:25:2925] % velocity range. The routine is parallelised around this velocity list

phi0 = -0.136 % initial latitude, this has very minimal effect
% zeta0 = 1.69 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
 zeta0 = 1.78 % initial heading angle, this is the zeta to reach close to 1.704 rad heading angle (SSO)
% for k = [33000:1000:36000] % altitude range
for k = [33000:500:36000]
% for j = [0 0.0125 0.025 0.0375 0.05] % trajectory angle range
% for j = [deg2rad(2.5): deg2rad(1): deg2rad(3.5)]
  for j = [deg2rad(1.5)]  
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

for i = 1:length(u)
% for i = 1:length(u)
u(i)

[mpayload(i), x, zeta, phi,Alt,v,t,Alpha,m,gamma,q,Vec_angle,T,CL,L,AltF(i),vF(i)] = ThirdStageOptm(k,j,u(i), phi0, zeta0, 0);

temp_payload(i) = mpayload(i);
end
mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',temp_payload.']]
temp_guess_no = temp_guess_no + 1;
end
end

dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')
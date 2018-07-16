function output = CombinedEndpoint(input)

t01 = input.phase(1).initialtime;
tf1 = input.phase(1).finaltime;
x01 = input.phase(1).initialstate;
xf1 = input.phase(1).finalstate;

t021 = input.phase(2).initialtime;
tf21 = input.phase(2).finaltime;
x021 = input.phase(2).initialstate;
xf21 = input.phase(2).finalstate;

t022 = input.phase(3).initialtime;
tf22 = input.phase(3).finaltime;
x022 = input.phase(3).initialstate;
xf22 = input.phase(3).finalstate;

%%

output.eventgroup(1).event = [x021(1)-xf1(1) x021(2)-xf1(9) x021(3)-xf1(8) x021(4)-xf1(2) x021(5)-xf1(4) x021(6)-xf1(6) x021(7)-xf1(5) t021-tf1];


%%

lon10  = input.phase(2).initialstate(2);
lat10  = input.phase(2).initialstate(3);
zeta10 = input.phase(2).initialstate(6);

% Calculate approximate initial launch site
initial_launch_lon = lon10-0.0051*cos(zeta10);
initial_launch_lat = lat10-0.0051*sin(zeta10);

alt1F  = input.phase(2).finalstate(1);
v1F    = input.phase(2).finalstate(4);
gamma1F  = input.phase(2).finalstate(5);

% const = input.auxdata.const;

% cost = input.phase(2).integral;

%%


lon22f  = input.phase(3).finalstate(2);
lat22f  = input.phase(3).finalstate(3);

output.eventgroup(2).event = [x022(1:9)-xf21(1:9)];
output.eventgroup(3).event = [t022-tf21];

%%
altf3 = input.phase(4).finalstate(1);
vf3 = input.phase(4).finalstate(2);
gammaf3 = input.phase(4).finalstate(3);
mf3 = input.phase(4).finalstate(4);
phif3 = input.phase(4).finalstate(6);
zetaf3 = input.phase(4).finalstate(7);
%%

t03 = input.phase(4).initialtime;
tf3 = input.phase(4).finaltime;
x03 = input.phase(4).initialstate;
xf3 = input.phase(4).finalstate;


[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff] = ThirdStageSim(altf3,gammaf3,vf3, phif3, 0,zetaf3, mf3, input.auxdata);

output.eventgroup(4).event =[xf21(1)-x03(1) xf21(4)-x03(2) xf21(5)-x03(3) xf21(6)-x03(7) xf21(3)-x03(6) t03-tf21 tf3-t03 AltF_actual inc_diff];



%% cost
output.objective = -mpayload;
% output.objective = cost;

end

function output = CombinedEndpoint(input)

V1f = input.phase(1).finalstate(1);
v1f = input.phase(1).finalstate(2); 
gamma1f = input.phase(1).finalstate(3); 
mfuel1f = input.phase(1).finalstate(4); 
gammadot1f = input.phase(1).finalstate(5);
zeta1f = input.phase(1).finalstate(6);
t1f= input.phase(1).finaltime;
%
t02= input.phase(2).initialtime;
alt02 = input.phase(2).initialstate(1);
v02 = input.phase(2).initialstate(2);
gamma02 = input.phase(2).initialstate(3);
m02 = input.phase(2).initialstate(4);

%

output.eventgroup(1).event = [alt02-V1f, v02-v1f, gamma02-gamma1f, t02-t1f];

%
altf2 = input.phase(2).finalstate(1);
vf2 = input.phase(2).finalstate(2);
gammaf2 = input.phase(2).finalstate(3);
mf2 = input.phase(2).finalstate(4);

%%

[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff]=ThirdStageSim(altf2,gammaf2,vf2, 0, 0,deg2rad(97), mf2, input.auxdata);

output.eventgroup(2).event = AltF_actual;

output.objective = -mpayload;

end

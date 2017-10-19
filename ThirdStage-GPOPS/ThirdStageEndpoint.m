function output = ThirdStageEndpoint(input)

altf = input.phase.finalstate(1);
vf = input.phase.finalstate(2);
gammaf = input.phase.finalstate(3);
mf = input.phase.finalstate(4);

%%

[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff]=ThirdStageSim(altf,gammaf,vf, 0, 0,deg2rad(97), mf, input.auxdata);

%% cost
output.objective = -mpayload;

output.eventgroup.event = AltF_actual;


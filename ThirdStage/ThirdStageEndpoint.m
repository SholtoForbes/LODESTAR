function output = ThirdStageEndpoint(input)

altf = input.phase.finalstate(1);
vf = input.phase.finalstate(2);
gammaf = input.phase.finalstate(3);
mf = input.phase.finalstate(4);
phif = input.phase.finalstate(6);
zetaf = input.phase.finalstate(7);
%%

[AltF_actual, vF, Alt, v, t, mpayload]=ThirdStageSim(altf,gammaf,vf, phif, 0,zetaf, mf, input.auxdata);

%% cost
output.objective = -mpayload;

output.eventgroup.event = AltF_actual;


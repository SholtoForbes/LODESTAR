function output = rlvEntryEndpoint(input)

% Inputs
% input.phase(phasenumber).initialstate -- row
% input.phase(phasenumber).finalstate -- row
% input.phase(phasenumber).initialtime -- scalar
% input.phase(phasenumber).finaltime -- scalar
% input.phase(phasenumber).integral -- row
%
% input.parameter -- row

% input.auxdata = auxiliary information

% Output
% output.objective -- scalar
% output.eventgroup(eventnumber).event -- row
% latf = input.phase.finalstate(3);

altf = input.phase.finalstate(1);

% vf = input.phase.finalstate(4);
vf = input.phase.finalstate(2);
gammaf = input.phase.finalstate(3);
mf = input.phase.finalstate(4);
% m    = input.phase.state(:,4);

%%

[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff]=ThirdStageSim(altf,gammaf,vf, 0, 0,deg2rad(97), mf, input.auxdata);

%% cost
% output.objective = -latf;
% output.objective = -vf;
% output.objective = (radf - input.auxdata.Re - 80000)^2;
%  output.objective = 0;
output.objective = -mpayload;
%  output.objective = -altf;

output.eventgroup.event = AltF_actual;


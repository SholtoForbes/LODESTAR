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
latf = input.phase.finalstate(3);
lonf = input.phase.finalstate(2);
vf = input.phase.finalstate(4);
mFuel0 = input.phase.initialstate(9);
% altf = input.phase.finalstate(2);
% cost
% output.objective = 1000*(latf+0.1)^2;
% output.objective = latf;
% output.objective = -vf;
output.objective = mFuel0;
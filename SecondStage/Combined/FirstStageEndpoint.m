function output = goddardRocketEndpoint(input)

auxdata = input.auxdata;

% Variables at Start and Terminus of Phase
t0 = input.phase.initialtime;
tf = input.phase.finaltime;
x0 = input.phase.initialstate;
xf = input.phase.finalstate;

% Objective Function
% output.objective = x0(3);
output.objective = -xf(2);

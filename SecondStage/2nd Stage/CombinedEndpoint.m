function output = CombinedEndpoint(input)

t021 = input.phase.initialtime;
tf21 = input.phase.finaltime;
x021 = input.phase.initialstate;
xf21 = input.phase.finalstate;


output.eventgroup.event = [tf21-t021];

cost = input.phase.integral;

%% cost
% output.objective = -mpayload;
output.objective = cost;

end

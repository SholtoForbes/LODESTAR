function output = chemicalReactionEndpoint(input)

output.objective = -input.phase.finalstate(2);

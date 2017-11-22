function output = SecondStageEndpoint(input)

alt1F  = input.phase(1).finalstate(1);
v1F    = input.phase(1).finalstate(4);
gamma1F  = input.phase(1).finalstate(5);

const = input.auxdata.const;

%% THIRD STAGE Payload Matrix %%===========================================
% NEED TO WATCH THIS, IT CAN EXTRAPOLATE BUT IT DOESNT DO IT WELL

if v1F > 2575
    ThirdStagePayloadMass = input.auxdata.PayloadGrid(alt1F,gamma1F,v1F);
else
    ThirdStagePayloadMass = gaussmf(v1F, [300 2575] )*input.auxdata.PayloadGrid(alt1F,gamma1F,2575);
end

output.objective = -ThirdStagePayloadMass;


%%
end

function output = CombinedEndpoint(input)

lon10  = input.phase(1).initialstate(2);
lat10  = input.phase(1).initialstate(3);
zeta10 = input.phase(1).initialstate(6);

% Calculate approximate initial launch site
initial_launch_lon = lon10-0.0051*cos(zeta10);
initial_launch_lat = lat10-0.0051*sin(zeta10);

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
t01 = input.phase(1).initialtime;
tf1 = input.phase(1).finaltime;
x01 = input.phase(1).initialstate;
xf1 = input.phase(1).finalstate;

t02 = input.phase(2).initialtime;
tf2 = input.phase(2).finaltime;
x02 = input.phase(2).initialstate;
xf2 = input.phase(2).finalstate;

lon2f  = input.phase(2).finalstate(2);
lat2f  = input.phase(2).finalstate(3);

output.eventgroup(1).event = [x02(1:9)-xf1(1:9) t02-tf1];
output.eventgroup(2).event = [lon2f-initial_launch_lon lat2f-initial_launch_lat];
end

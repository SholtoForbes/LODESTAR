function output = FirstStageContinuous(input)


auxdata = input.auxdata;

t = input.phase.time';
x = input.phase.state';
u = input.phase.control';

phase = 'postpitch';

interp = auxdata.interp;
Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

[dz,q,xi] = rocketDynamics(x,u,t,phase,interp,Throttle,Vehicle,Atmosphere); % Pass primal variables into dynamics simulator


output.path = q';


output.dynamics = dz';

%------------------------------------%
% End File:  goddardRocketContinuous %
%------------------------------------%

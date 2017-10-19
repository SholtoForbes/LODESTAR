function phaseout = ThirdStageContinuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%

alt  = input.phase.state(:,1);
v    = input.phase.state(:,2);
gamma  = input.phase.state(:,3);
m    = input.phase.state(:,4);
Alpha    = input.phase.state(:,5);
time = input.phase.time;

Alphadot  = input.phase.control(:,1);

auxdata=input.auxdata;

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T] = ThirdStageDyn(alt,gamma,v,m,Alpha,time,auxdata, Alphadot);

phaseout.dynamics  = [rdot.', vdot.', gammadot.', -mdot*ones(length(rdot),1), Alphadot];

Alpha_constraint = Alpha-AoA_max;

phaseout.path = [Vec_angle,Alpha_constraint];
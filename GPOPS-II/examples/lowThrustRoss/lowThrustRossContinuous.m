function phaseout = lowThrustRossContinuous(input)

T = input.auxdata.T;
m0 = input.auxdata.m0;
dm = input.auxdata.dm;
mu = input.auxdata.mu;

t      = input.phase.time;
r      = input.phase.state(:,1);
theta  = input.phase.state(:,2);
vr     = input.phase.state(:,3);
vtheta = input.phase.state(:,4);
ur     = input.phase.control(:,1);
ut     = input.phase.control(:,2);

dr      = vr;
dtheta  = vtheta./r;
dvr     = (vtheta.^2)./r-mu./(r.^2)+ur;
dvtheta = -(vr.*vtheta)./r+ut;

phaseout.dynamics  = [dr, dtheta, dvr, dvtheta];

function phaseout = chemicalReactionContinuous(input)

rho = input.auxdata.rho;
k   = input.auxdata.k;

t  = input.phase.time;
x1 = input.phase.state(:,1);
x2 = input.phase.state(:,2);
u  = input.phase.control(:,1);

x1dot = -u.*x1;
x2dot = -x1dot - rho.*(u.^k).*x2;

phaseout.dynamics = [x1dot, x2dot];

function phaseout = ThirdStageContinuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
% r  = input.phase.state(:,1);
% xi  = input.phase.state(:,2);
% phi  = input.phase.state(:,3);
% v    = input.phase.state(:,4);
% gamma  = input.phase.state(:,5);
% zeta  = input.phase.state(:,6);
% m    = input.phase.state(:,7);


alt  = input.phase.state(:,1);
v    = input.phase.state(:,2);
gamma  = input.phase.state(:,3);
m    = input.phase.state(:,4);
Alpha    = input.phase.state(:,5);

xi = 0;
phi = 0;
zeta = deg2rad(97);

Alphadot  = input.phase.control(:,1);
% bank = input.phase.control(:,2);

% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%
% cd0      = input.auxdata.cd(1);
% cd1      = input.auxdata.cd(2);
% cd2      = input.auxdata.cd(3);
% cl0      = input.auxdata.cl(1);
% cl1      = input.auxdata.cl(2);
% mu       = input.auxdata.mu;
% rho0     = input.auxdata.rho0;
% H        = input.auxdata.H;
% S        = input.auxdata.S;
% mass     = input.auxdata.mass;
% altitude = rad - input.auxdata.Re;
% CD       = cd0+cd1*aoa+cd2*aoa.^2;
% rho      = rho0*exp(-altitude/H);
% CL       = cl0+cl1*aoa;
% q        = 0.5*rho.*v.^2;
% D        = q.*S.*CD./mass;
% L        = q.*S.*CL./mass;
% gravity  = mu./rad.^2;
% 
% % ---------------------------------------------------%
% % ---- Evaluate Right-Hand Side of the Dynamics ---- %
% % ---------------------------------------------------%
% raddot = v.*sin(fpa);
% londot = v.*cos(fpa).*sin(azi)./(rad.*cos(lat));
% latdot = v.*cos(fpa).*cos(azi)./rad;
% vdot   = -D-gravity.*sin(fpa);
% fpadot = (L.*cos(bank)-cos(fpa).*(gravity-v.^2./rad))./v;
% azidot = (L.*sin(bank)./cos(fpa)+v.^2.*cos(fpa).*sin(azi).*tan(lat)./rad)./v;

auxdata=input.auxdata;

[rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T] = ThirdStageDyn(alt,xi,phi,gamma,v,zeta,m,Alpha,auxdata);

% Vec_angle
% phaseout.dynamics  = [rdot, xidot, phidot, vdot, gammadot, zetadot, -mdot*ones(length(rdot),1)];
phaseout.dynamics  = [rdot, vdot, gammadot, -mdot*ones(length(rdot),1), Alphadot];
% AoA_max
% Alpha
Alpha_constraint = Alpha-AoA_max;
% Vec_angle
phaseout.path = [Vec_angle,Alpha_constraint];
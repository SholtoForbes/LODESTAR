function phaseout = rlvEntryContinuous(input)
% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
global rad
rad  = input.phase.state(:,1);
lon  = input.phase.state(:,2);
lat  = input.phase.state(:,3);
v    = input.phase.state(:,4);
fpa  = input.phase.state(:,5);
azi  = input.phase.state(:,6);

aoa  = input.phase.state(:,7);
bank  = input.phase.state(:,8);

% aoa  = input.phase.control(:,1);
% bank = input.phase.control(:,2);

aoadot  = input.phase.control(:,1);
bankdot = input.phase.control(:,2);

% bank = 0*ones(length(aoa),1);
% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

global time
time = input.phase(1).time;

auxdata = input.auxdata;

[raddot,londot,latdot,fpadot,vdot,azidot, q, M, Fd, rho,L] = VehicleModelReturn(fpa, rad, v,auxdata,azi,lat,lon,aoa,bank);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout.dynamics  = [raddot.', londot.', latdot.', vdot.', fpadot.', azidot.', aoadot, bankdot];
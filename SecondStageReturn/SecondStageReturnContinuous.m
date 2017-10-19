function phaseout = SecondStageReturnContinuous(input)
% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rad  = input.phase.state(:,1);
lon  = input.phase.state(:,2);
lat  = input.phase.state(:,3);
v    = input.phase.state(:,4);
fpa  = input.phase.state(:,5);
azi  = input.phase.state(:,6);

aoa  = input.phase.state(:,7);
bank  = input.phase.state(:,8);
mFuel  = input.phase.state(:,9);
throttle  = input.phase.state(:,10);

% aoa  = input.phase.control(:,1);
% bank = input.phase.control(:,2);

aoadot  = input.phase.control(:,1);
bankdot = input.phase.control(:,2);
throttledot = input.phase.control(:,3);

% bank = 0*ones(length(aoa),1);
% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time = input.phase(1).time;

auxdata = input.auxdata;

[raddot,londot,latdot,fpadot,vdot,azidot, q, M, Fd, rho,L,Fueldt,T] = VehicleModelReturn(fpa, rad, v,auxdata,azi,lat,lon,aoa,bank,throttle, mFuel);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

% phaseout.dynamics  = [raddot.', londot.', latdot.', vdot.', fpadot.', azidot.', aoadot, bankdot];
phaseout.dynamics  = [raddot, londot, latdot, vdot, fpadot, azidot, aoadot, bankdot, -Fueldt, throttledot];

% Thrust_Constraint = Fueldt;
% Thrust_Constraint(M>5.1) = 0;

% Throttle_Constraint = throttle;
% Throttle_Constraint(M>5.1) = 0;

phaseout.path = q;
% phaseout.path = [q, Throttle_Constraint];
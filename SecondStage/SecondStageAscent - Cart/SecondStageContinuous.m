function phaseout = SecondStageContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics
alt1  = input.phase(1).state(:,1);
lon1  = input.phase(1).state(:,2);
lat1  = input.phase(1).state(:,3);
v1    = input.phase(1).state(:,4);
gamma1  = input.phase(1).state(:,5);
zeta1  = input.phase(1).state(:,6);

Alpha1  = input.phase(1).state(:,7);
% eta1  = input.phase(1).state(:,8);
mFuel1  = input.phase(1).state(:,8);
throttle1  = 1;
eta1 = 0;

aoadot1  = input.phase(1).control(:,1);
% bankdot1 = input.phase(1).control(:,2);

% bank = 0*ones(length(aoa),1);
% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time1 = input.phase(1).time;

auxdata = input.auxdata;

[altdot1,londot1,latdot1,fpadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T] = VehicleModelAscent(gamma1, alt1, v1,auxdata,zeta1,lat1,lon1,Alpha1,eta1,throttle1, mFuel1,mFuel1(1),mFuel1(end), 1);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(1).dynamics  = [altdot1, londot1, latdot1, vdot1, fpadot1, azidot1, aoadot1, -Fueldt1];
phaseout(1).path = q1;


end

%======================================================
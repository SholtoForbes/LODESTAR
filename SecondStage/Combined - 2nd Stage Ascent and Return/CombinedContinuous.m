function phaseout = CombinedContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics
alt1  = input.phase(1).state(:,1);
lon1  = input.phase(1).state(:,2);
lat1  = input.phase(1).state(:,3);
v1    = input.phase(1).state(:,4);
gamma1  = input.phase(1).state(:,5);
zeta1  = input.phase(1).state(:,6);

Alpha1  = input.phase(1).state(:,7);
eta1  = input.phase(1).state(:,8);
mFuel1  = input.phase(1).state(:,9);
throttle1  = 1;

aoadot1  = input.phase(1).control(:,1);
bankdot1 = input.phase(1).control(:,2);

% bank = 0*ones(length(aoa),1);
% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time1 = input.phase(1).time;

auxdata = input.auxdata;

[altdot1,londot1,latdot1,fpadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T] = VehicleModelCombined(gamma1, alt1, v1,auxdata,zeta1,lat1,lon1,Alpha1,eta1,throttle1, mFuel1, 1);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(1).dynamics  = [altdot1, londot1, latdot1, vdot1, fpadot1, azidot1, aoadot1, bankdot1, -Fueldt1];
phaseout(1).path = q1;
%%

alt2  = input.phase(2).state(:,1);
lon2  = input.phase(2).state(:,2);
lat2  = input.phase(2).state(:,3);
v2    = input.phase(2).state(:,4);
fpa2  = input.phase(2).state(:,5);
azi2  = input.phase(2).state(:,6);

aoa2  = input.phase(2).state(:,7);
bank2  = input.phase(2).state(:,8);
mFuel2  = input.phase(2).state(:,9);
throttle2  = input.phase(2).state(:,10);

aoadot2  = input.phase(2).control(:,1);
bankdot2 = input.phase(2).control(:,2);
throttledot2 = input.phase(2).control(:,3);

% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time2 = input.phase(1).time;

[altdot2,londot2,latdot2,fpadot2,vdot2,azidot2, q2, M, Fd, rho,L,Fueldt2,T] = VehicleModelCombined(fpa2, alt2, v2,auxdata,azi2,lat2,lon2,aoa2,bank2,throttle2, mFuel2,0);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(2).dynamics  = [altdot2, londot2, latdot2, vdot2, fpadot2, azidot2, aoadot2, bankdot2, -Fueldt2, throttledot2];


phaseout(2).path = q2;

end

%======================================================
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

[altdot1,londot1,latdot1,fpadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T3,Isp1,q1_aftershock] = VehicleModelCombined(gamma1, alt1, v1,auxdata,zeta1,lat1,lon1,Alpha1,eta1,throttle1, mFuel1,mFuel1(1),mFuel1(end), 1);

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

[altdot2,londot2,latdot2,fpadot2,vdot2,azidot2, q2, M, Fd, rho,L,Fueldt2,T3,Isp2,q2_aftershock] = VehicleModelCombined(fpa2, alt2, v2,auxdata,azi2,lat2,lon2,aoa2,bank2,throttle2, mFuel2,0,0,0);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(2).dynamics  = [altdot2, londot2, latdot2, vdot2, fpadot2, azidot2, aoadot2, bankdot2, -Fueldt2, throttledot2];


phaseout(2).path = q2;

%% Third Stage

alt3  = input.phase(3).state(:,1);
v3    = input.phase(3).state(:,2);
gamma3  = input.phase(3).state(:,3);
m3    = input.phase(3).state(:,4);
Alpha3    = input.phase(3).state(:,5);
phi3    = input.phase(3).state(:,6);
zeta3    = input.phase(3).state(:,7);
time3 = input.phase(3).time;

Alphadot  = input.phase(3).control(:,1);

auxdata=input.auxdata;

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3] = ThirdStageDyn(alt3,gamma3,v3,m3,Alpha3,time3,auxdata, Alphadot, phi3, zeta3);

phaseout(3).dynamics  = [rdot3, vdot3, gammadot3, -mdot3*ones(length(rdot3),1), Alphadot, phidot3, zetadot3];

Alpha_constraint = Alpha3-AoA_max3;

phaseout(3).path = [Vec_angle3,Alpha_constraint];


end

%======================================================
function phaseout = CombinedContinuous(input)
%% First Stage

auxdata = input.auxdata;

t = input.phase(1).time';
x = input.phase(1).state';
u = input.phase(1).control';

phase = 'postpitch';

interp = auxdata.interp;
Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

[dz1,q1,xi1] = rocketDynamics(x,u,t,phase,interp,Throttle,Vehicle,Atmosphere); % Pass primal variables into dynamics simulator


phaseout(1).path = q1';


phaseout(1).dynamics = dz1';



%% Second Stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alt21  = input.phase(2).state(:,1);
lon21  = input.phase(2).state(:,2);
lat21  = input.phase(2).state(:,3);
v21    = input.phase(2).state(:,4);
gamma21  = input.phase(2).state(:,5);
zeta21  = input.phase(2).state(:,6);

Alpha21  = input.phase(2).state(:,7);
eta21  = input.phase(2).state(:,8);
mFuel21  = input.phase(2).state(:,9);
throttle21  = 1;

aoadot21  = input.phase(2).control(:,1);
bankdot21 = input.phase(2).control(:,2);

% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time21 = input.phase(2).time;



[altdot21,londot21,latdot21,fpadot21,vdot21,azidot21, q21, M, Fd, rho,L,Fueldt21,T3,Isp21,q21_aftershock] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,Alpha21,eta21,throttle21, mFuel21,mFuel21(1),mFuel21(end), 1, 0);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(2).dynamics  = [altdot21, londot21, latdot21, vdot21, fpadot21, azidot21, aoadot21, bankdot21, -Fueldt21];
phaseout(2).path = q21;
%%

alt22  = input.phase(3).state(:,1);
lon22  = input.phase(3).state(:,2);
lat22  = input.phase(3).state(:,3);
v22    = input.phase(3).state(:,4);
fpa22  = input.phase(3).state(:,5);
azi22  = input.phase(3).state(:,6);

aoa22  = input.phase(3).state(:,7);
bank22  = input.phase(3).state(:,8);
mFuel22  = input.phase(3).state(:,9);
throttle22  = input.phase(3).state(:,10);

aoadot22  = input.phase(3).control(:,1);
bankdot22 = input.phase(3).control(:,2);
throttledot22 = input.phase(3).control(:,3);

% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time22 = input.phase(2).time;

[altdot22,londot22,latdot22,fpadot22,vdot22,azidot22, q22, M, Fd, rho,L,Fueldt22,T22,Isp22,q22_aftershock] = VehicleModelCombined(fpa22, alt22, v22,auxdata,azi22,lat22,lon22,aoa22,bank22,throttle22, mFuel22,0,0,0, 0);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout(3).dynamics  = [altdot22, londot22, latdot22, vdot22, fpadot22, azidot22, aoadot22, bankdot22, -Fueldt22, throttledot22];


phaseout(3).path = q22;

%% Third Stage

alt3  = input.phase(4).state(:,1);
v3    = input.phase(4).state(:,2);
gamma3  = input.phase(4).state(:,3);
m3    = input.phase(4).state(:,4);
Alpha3    = input.phase(4).state(:,5);
phi3    = input.phase(4).state(:,6);
zeta3    = input.phase(4).state(:,7);
time3 = input.phase(4).time;

Alphadot  = input.phase(4).control(:,1);

auxdata=input.auxdata;

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3] = ThirdStageDyn(alt3,gamma3,v3,m3,Alpha3,time3,auxdata, Alphadot, phi3, zeta3);

phaseout(4).dynamics  = [rdot3, vdot3, gammadot3, -mdot3*ones(length(rdot3),1), Alphadot, phidot3, zetadot3];

Alpha_constraint = Alpha3-AoA_max3;

phaseout(4).path = [Vec_angle3,Alpha_constraint];


end

%======================================================
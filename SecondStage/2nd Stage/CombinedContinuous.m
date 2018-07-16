function phaseout = CombinedContinuous(input)


%% Second Stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alt21  = input.phase.state(:,1);
lon21  = input.phase.state(:,2);
lat21  = input.phase.state(:,3);
v21    = input.phase.state(:,4);
gamma21  = input.phase.state(:,5);
zeta21  = input.phase.state(:,6);

Alpha21  = input.phase.state(:,7);
eta21  = input.phase.state(:,8);
mFuel21  = input.phase.state(:,9);
throttle21  = 1;

aoadot21  = input.phase.control(:,1);
bankdot21 = input.phase.control(:,2);
% throttle21 = input.phase.control(:,3);
% ---------------------------------------------------%
% ------- Compute the Aerodynamic Quantities --------%
% ---------------------------------------------------%

time21 = input.phase.time;

auxdata = input.auxdata;

[altdot21,londot21,latdot21,fpadot21,vdot21,azidot21, q21, M, Fd, rho,L,Fueldt21,T3,Isp21,q21_aftershock] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,Alpha21,eta21,throttle21, mFuel21,mFuel21(1),mFuel21(end), 1, 0);

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

phaseout.dynamics  = [altdot21, londot21, latdot21, vdot21, fpadot21, azidot21, aoadot21, bankdot21, -Fueldt21];
phaseout.path = q21;

phaseout.integrand = (q21-50000).^2;


end

%======================================================
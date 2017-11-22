function phaseout = SecondStageContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alt = input.phase.state(:,1).';
v = input.phase.state(:,2).'; 
gamma = input.phase.state(:,3).'; 
mfuel = input.phase.state(:,4).'; 
zeta = input.phase.state(:,5).';
Alpha = input.phase.state(:,6).';

Alphadot  = input.phase.control; 

time = input.phase.time.';

Stage2 = input.auxdata.Stage2;
Stage3 = input.auxdata.Stage3;
interp = input.auxdata.interp;
const = input.auxdata.const;

auxdata = input.auxdata;

[altdot, dfuel, Fueldt, vdot, q, M, D, Thrust, flapdeflection, rho,zeta,phi,eq,zetadot,xi, gammadot] = VehicleModel(time, gamma, alt, v, mfuel,interp,const, interp.Atmosphere,zeta,auxdata,Alpha,gamma);


mfueldot = -Fueldt;
%==========================================================================
% Thrust
% D
% M
% vdot.'
phaseout.dynamics = [altdot.',vdot.', gammadot.', mfueldot.', zetadot.',Alphadot];

phaseout.path = [q.'];
end

%======================================================
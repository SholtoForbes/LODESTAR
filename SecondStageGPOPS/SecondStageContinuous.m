function phaseout = SecondStageContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics

% This uses velocity calculated in the Cost file

% This file is calculated after the Cost file in the iterative process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = input.phase.state(:,1).';
v = input.phase.state(:,2).'; 
gamma = input.phase.state(:,3).'; 
mfuel = input.phase.state(:,4).'; 
gammadot = input.phase.state(:,5).';
zeta = input.phase.state(:,6).';

omegadot  = input.phase.control.'; 

time = input.phase.time.';

Stage2 = input.auxdata.Stage2;
Stage3 = input.auxdata.Stage3;
interp = input.auxdata.interp;
const = input.auxdata.const;

auxdata = input.auxdata;

[dfuel, Engine.Fueldt, a, q, M, Vehicle.Fd, Engine.Thrust, Vehicle.flapdeflection, Vehicle.Alpha, rho,Vehicle.lift,zeta,phi,Engine.eq,zetadot] = VehicleModel(time, gamma, V, v, mfuel,interp,const,gammadot, interp.Atmosphere,zeta,Stage2.mStruct,Stage3.mTot,auxdata);

vdot = a;
mfueldot = -Engine.Fueldt; 
%==========================================================================

Vdot = v.*sin(gamma);

%==========================================================================

phaseout.dynamics = [Vdot.',vdot.', gammadot.', mfueldot.', omegadot.', zetadot.'];

phaseout.path = [q.',Vehicle.Alpha.'];
end

%======================================================
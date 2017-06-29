%% Program to simulate a banking turn
% Sholto Forbes-Spyratos
clear all
mSPARTAN_empty = 4910.5 - 132.8 + 179.41; % Mass of the empty SPARTAN, with scaled fuel tank mass

V0 = 36800;
gamma0 = 0.048;
zeta0 = 1.69;
phi0 = -0.12913;
xi0 = 0;
v0 = 2920;

Initial_States = [V0,phi0,gamma0,v0,zeta0];


communicator = importdata('communicator.txt');
communicator_trim = importdata('communicator_trim.txt');
Atmosphere = dlmread('atmosphere.txt');

interp.flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
interp.flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,3),communicator_trim(:,5));
interp.flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,3),communicator_trim(:,6));

[MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';

interp.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
interp.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
interp.pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');

Alpha = 7; % aoa (deg)
FlapDeflection = 0;
eta = 1; % roll (rad)

[t, y] = ode45(@(f_t,f_y) ForwardSimReturn(f_y,Alpha,eta,Atmosphere,interp,FlapDeflection,mSPARTAN_empty),[0 400],Initial_States);


figure(401)
subplot(5,5,[1,5]);
plot(t,y(:,1));
title('Altitude');

subplot(5,5,[6,10]);
plot(t,y(:,5));
title('Heading Angle');

subplot(5,5,[11,15]);
plot(t,y(:,4));
title('Velocity');

subplot(5,5,[16,20]);
plot(t,y(:,2));
title('Latitude');











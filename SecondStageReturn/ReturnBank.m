1%% Program to simulate a banking turn
% Sholto Forbes-Spyratos
clear all
mSPARTAN_empty = 4910.5 - 132.8 + 179.41; % Mass of the empty SPARTAN, with scaled fuel tank mass

V0 = 35000;
gamma0 = 0.048;
zeta0 = 1.76;
phi0 = -0.12913;
xi0 = 0;
v0 = 2920;

Initial_States = [V0,phi0,gamma0,v0,zeta0,xi0];

%% Define Latitude and Longitude Constraints
% Define the latitude and longitude of home base
returncond.phi = -0.15;
% returncond.xi = 0;

% Define acceptable landing radius
returncond.rad = 0.05; %in radians

%%

% communicator = importdata('communicator.txt');
% communicator_trim = importdata('communicator_trim.txt');

aero = importdata('SPARTANaero.txt');

Atmosphere = dlmread('atmosphere.txt');

% interp.flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
% interp.flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,3),communicator_trim(:,5));
% interp.flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,3),communicator_trim(:,6));

% [MList,AOAList] = ndgrid(unique(aero(:,1)),unique(aero(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% pitchingmoment_Grid = reshape(aero(:,11),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% 
% interp.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
% interp.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
% interp.pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');


% interp.Cl_spline = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,3));
% interp.Cd_spline = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,4));

aero = importdata('SPARTANaero.txt');

interp.Cl_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,3));
interp.Cd_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,4));


[MList,AOAList] = ndgrid(unique(aero(:,1)),unique(aero(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid = [];
Cd_Grid = [];

for i = 1:numel(MList)
    M_temp = MList(i);
    AoA_temp = AOAList(i);
    
    Cl_temp = interp.Cl_scattered(M_temp,AoA_temp);
    Cd_temp = interp.Cd_scattered(M_temp,AoA_temp);
    
    I = cell(1, ndims(MList)); 
    [I{:}] = ind2sub(size(MList),i);
    
    Cl_Grid(I{(1)},I{(2)}) = Cl_temp;
    Cd_Grid(I{(1)},I{(2)}) = Cd_temp;

end

interp.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
interp.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');


% Alpha = 7; % aoa (deg)
% FlapDeflection = 0;
% eta = 1; % roll (rad)

% [t, y] = ode45(@(f_t,f_y) ForwardSimReturn(f_y,Alpha,eta,Atmosphere,interp,FlapDeflection,mSPARTAN_empty),[0 400],Initial_States);

% options.Algorithm = 'sqp';
options.Display = 'iter';
% options.MaxFunEvals = 5000;

% options.ScaleProblem = 'obj-and-constr';
% options.DiffMinChange = 0.0005;

num_div = 20;% no of timestep divisions

Altitude_0 = V0-V0*(1:(num_div-1))/(num_div-1);

% controls0 = [7*ones(1,num_div) 0.5*ones(1,num_div) 450];
% lb = [1*ones(1,num_div) -1*ones(1,num_div) 400];
% ub = [8*ones(1,num_div) deg2rad(80)*ones(1,num_div) 600];


eta_00 = 1; % initial roll guess
AoA_00 = 7;

controls0 = [0*ones(1,num_div) 0*ones(1,num_div) 600 eta_00 AoA_00];

lb = [-.1*ones(1,num_div) -.01*ones(1,num_div) 400 -1 0];
ub = [.1*ones(1,num_div) .01*ones(1,num_div) 600 1.5 8];

[controls_opt,fval,exitflag] = fmincon(@(controls)BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,num_div),controls0,[],[],[],[],lb,ub,@(controls)ReturnConstraint(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,returncond,num_div),options)
% [controls_opt,fval,exitflag] = particleswarm(@(controls)BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,num_div),length(lb),lb,ub)

[cost,phi,t,y,q,xi,zeta] = BankOpt(controls_opt,Initial_States,Atmosphere,interp,mSPARTAN_empty,num_div,options);


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

% subplot(5,7,[21,25]);
% plot(t,y(:,2));
% title('Latitude');
% 
% subplot(5,7,[26,30]);
% plot(t,y(:,2));
% title('Latitude');











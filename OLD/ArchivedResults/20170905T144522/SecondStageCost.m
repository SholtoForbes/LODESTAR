function [EndpointCost, RunningCost] = TwoStage2DCost(primal, algorithm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function for Rocket-Scramjet-Rocket System

% This module contains the entire vehicle/system model, as well as defining
% the final cost

% The previous iteration of velocity is initially calculated using primals, and subsequently
% dynamic calculations are performed to produce next velocity step. These
% will eventually converge.

% Global variables are used to pass this velocity to the Dynamics file


%% Import Data %%==========================================================
% Import Primals
V = primal.states(1, :); % Scaled vertical position
v = primal.states(2,:);
gamma  = primal.states(3, :); % Velocity angle
mfuel = primal.states(4,:);
gammadot = primal.states(5,:);
zeta = primal.states(6,:);
time = primal.nodes(1, :);

% Import Globals
global nodes
%% Vehicle Model %%========================================================
% This is the main dynamic simulation of the vehicle.
global M
global q
global dfuel
global a
global Fd
global Fueldt
global Thrust
global lift
global flapdeflection
global Alpha
global interp
global rho
global const
global gridded
global Atmosphere
global iteration
% global zeta
global phi
global eq
global zetadot
global mstruct
global mThirdStage
global Stage2
global Stage3

[dfuel, Fueldt, a, q, M, Fd, Thrust, flapdeflection, Alpha, rho,lift,zeta,phi,eq,zetadot] = VehicleModel(time, gamma, V, v, mfuel, nodes,interp,gridded,const,gammadot, interp.Atmosphere,zeta,Stage2.mStruct,Stage3.mTot);

IspNet = (Thrust-Fd)./Fueldt./9.81;
%% THIRD STAGE Payload Matrix %%===========================================
% NEED TO WATCH THIS, IT CAN EXTRAPOLATE BUT IT DOESNT DO IT WELL

global PayloadGrid
if v(end) > 2850
    ThirdStagePayloadMass = PayloadGrid(V(end),gamma(end),v(end));
else
    ThirdStagePayloadMass = gaussmf(v(end), [300 2850] )*PayloadGrid(V(end),gamma(end),2850);
end

% Define Cost =======================================================
% The pseudospectral solver is able to run with no cost at all, this is useful for checking dynamics

if const == 3  || const == 32
Endcost = 0;
% Endcost = 0.001*gamma(end);
% Endcost =1*sum(abs(omegadot));
elseif const == 1  || const == 12 || const == 13 || const == 14  || const == 15
Endcost =  - 0.01*mfuel(end) - ThirdStagePayloadMass;
end

EndpointCost = Endcost;


if const == 3 
    RunningCost =((q-50000).^2+10000)/10000 ; %
elseif const == 32
    RunningCost =((q-50000).^2+10000)/10000;
elseif const == 1  || const == 12 || const == 13 || const == 14  || const == 15
    RunningCost = 0;
end


%% Display Progress %%=====================================================
iteration = iteration + 1;
global iterative_V
global iterative_t
global iterative_V_f
if rem(iteration,5000) == 0
    iterative_V_f(end+1,:) = cumtrapz(time,v.*sin(gamma));
    cla
    hold on
    ThirdStagePayloadMass
    todisp = ['Total Iterations Done:',num2str(iteration)];
    disp(todisp);
    figure(10)
    iterative_V(end+1,:) = V;
    iterative_t(end+1,:) = time;
    if length(iterative_V(:,1)) > 4
        plot(iterative_t(end-4,:),iterative_V(end-4,:),'Color',[0.8 0.8 0.8])
    end
    if length(iterative_V(:,1)) > 2
        plot(iterative_t(end-2,:),iterative_V(end-2,:),'Color',[0.4 0.4 0.4])
    end
    if length(iterative_V(:,1)) > 3
        plot(iterative_t(end-3,:),iterative_V(end-3,:),'Color',[0.6 0.6 0.6])
    end
    plot(iterative_t(end-1,:),iterative_V(end-1,:),'Color',[0.2 0.2 0.2])
    plot(iterative_t(end,:),iterative_V(end,:),'Color',[0 0 0])
    filename = 'testnew51.gif';
    frame = getframe(10);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
    drawnow
end

end

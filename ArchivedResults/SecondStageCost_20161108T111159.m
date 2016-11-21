function [EndpointCost, RunningCost] = TwoStage2DCost(primal, algorithm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function for Rocket-Scramjet-Rocket System

% This module contains the entire vehicle/system model, as well as defining
% the final cost

% The previous iteration of velocity is initially calculated using primals, and subsequently
% dynamic calculations are performed to produce next velocity step. These
% will eventually converge.

% Global variables are used to pass this velocity to the Dynamics file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global nodes
% =======================================================
% Vehicle Model:
% =======================================================
% =========================================================================================
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
global scattered
global rho
global const
global scale
global grid
global Atmosphere
global iteration
global zeta
global phi

iteration = iteration + 1;

V = primal.states(1, :)*scale.V ; % Scaled vertical position
v = primal.states(2,:)*scale.v ;
theta  = primal.states(3, :)*scale.theta; % Velocity angle
mfuel = primal.states(4,:)*scale.m ;
time = primal.nodes(1, :); % Time
% thetadot  = primal.controls(1, :)*scale.theta;
thetadot = primal.states(5,:)*scale.thetadot;
omegadot  = primal.controls(1,:)*scale.thetadot;

[dfuel, Fueldt, a, q, M, Fd, Thrust, flapdeflection, Alpha, rho,lift,Penalty,zeta,phi] = VehicleModel(time, theta, V, v, mfuel, nodes,scattered,grid,const,thetadot, Atmosphere);

% THIRD STAGE ======================================================
% NEED TO WATCH THIS, IT CAN EXTRAPOLATE BUT IT DOESNT DO IT WELL

global ThirdStagePayloadMass
global alt_list
global gamma_list
global v_list
global payload_array
t1 = cputime;
if V(end) > 40000
    if theta(end) < 0
    ThirdStagePayloadMass = gaussmf(theta(end),[0.1 0])*gaussmf(V(end),[10000 40000])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,40000, 0, v(end));    
    else
    ThirdStagePayloadMass = gaussmf(V(end),[10000 40000])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,40000, theta(end), v(end));    
    end
elseif v(end) < 2500
    if theta(end) < 0
    ThirdStagePayloadMass = gaussmf(theta(end),[0.1 0])*gaussmf(v(end),[1000 2000])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), 2000);  
    else
    ThirdStagePayloadMass = gaussmf(v(end),[1000 2000])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), 2000);
    end
elseif V(end) < 34000
    if theta(end) < 0
    ThirdStagePayloadMass = gaussmf(theta(end),[0.1 0])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), v(end));  
    else
    ThirdStagePayloadMass = interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), v(end));
    end
elseif V(end) > 36500 && V(end) <= 40000
    if theta(end) < 0
    ThirdStagePayloadMass = gaussmf(theta(end),[0.1 0])*interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), v(end));  
    else
    ThirdStagePayloadMass = interp3(alt_list.course,gamma_list.course,v_list.course,payload_array.course,V(end), theta(end), v(end));
    end
elseif V(end) >=34000 && V(end) <= 36500
    if theta(end) < 0
    ThirdStagePayloadMass = gaussmf(theta(end),[0.1 0])*interp3(alt_list.fine,gamma_list.fine,v_list.fine,payload_array.fine,V(end), theta(end), v(end),'spline');  
    else
    ThirdStagePayloadMass = interp3(alt_list.fine,gamma_list.fine,v_list.fine,payload_array.fine,V(end), theta(end), v(end),'spline');
    end
end

time = cputime - t1;
% Define Cost =======================================================
% The pseudospectral solver is able to run with no cost at all, this is useful for checking dynamics

if const == 3 
Endcost = 0;
end

if const == 1  || const == 12 || const == 13 || const == 14
Endcost =  - 0.01*mfuel(end) - ThirdStagePayloadMass;
end

EndpointCost = Endcost;

if const == 1  || const == 12 || const == 13 || const == 14
    
%smoothing functions (can be adjusted depending on needs, remove if not working
% omegadot = diff(thetadot)./diff(time);
%  RunningCost = [0 0.005*abs(omegadot)]; % for smoothing 50kPa
% RunningCost = [0 0.001*abs(omegadot)]; %for smoothing 45kPa and 55kPa and
% high drag

    RunningCost = Penalty + abs(omegadot);
%     RunningCost = Penalty ; % The Penalty function ensures that it does not go over 50kPa, but still allows it to search that space. 
    %Omegadot cost smooths the trajectory (sometimes), but also throws a SOL error and
    %reduces the validity of the optimisation result (so only use if you
    %are sure of the end trajectory shape)
% RunningCost = 0;
end

if const == 3 

% RunningCost =((q-50000).^2+1000000)/1000000 + [0 0.005*abs(omegadot)]; % if a cost does not work, try loosening it 
% RunningCost =((q-50000).^2+2000000)/2000000 + [0 0.005*abs(omegadot)];
% RunningCost =((q-50000).^2+2000000)/2000000 + abs(thetadot);
% RunningCost =((q-50000).^2+2000000)/2000000; %looser
RunningCost =((q-50000).^2+10000)/10000; % tighter
end


end

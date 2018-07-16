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
global gridded
global Atmosphere
global iteration
% global zeta
global phi
global SPARTAN_SCALE
global eq
global zetadot
iteration = iteration + 1;


V = primal.states(1, :)*scale.V ; % Scaled vertical position
v = primal.states(2,:)*scale.v ;
theta  = primal.states(3, :)*scale.theta; % Velocity angle
mfuel = primal.states(4,:)*scale.m ;
time = primal.nodes(1, :); % Time
% thetadot  = primal.controls(1, :)*scale.theta;
thetadot = primal.states(5,:)*scale.thetadot;
zeta = primal.states(6,:);
omegadot  = primal.controls(1,:)*scale.thetadot;



[dfuel, Fueldt, a, q, M, Fd, Thrust, flapdeflection, Alpha, rho,lift,Penalty,zeta,phi,eq,zetadot] = VehicleModel(time, theta, V, v, mfuel, nodes,scattered,gridded,const,thetadot, Atmosphere, SPARTAN_SCALE,zeta);

IspNet = (Thrust-Fd)./Fueldt./9.81;
% THIRD STAGE ======================================================
% NEED TO WATCH THIS, IT CAN EXTRAPOLATE BUT IT DOESNT DO IT WELL


global PayloadGrid
% ThirdStagePayloadMass = PayloadGrid(phi(end),zeta(end),V(end),theta(end),v(end));
if v(end) > 2850
ThirdStagePayloadMass = PayloadGrid(V(end),theta(end),v(end));
else
    ThirdStagePayloadMass = gaussmf(v(end), [300 2850] )*PayloadGrid(V(end),theta(end),2850);
end

% Define Cost =======================================================
% The pseudospectral solver is able to run with no cost at all, this is useful for checking dynamics

if const == 3 
Endcost = 0;
Endcost =10*sum(abs(omegadot));
elseif const == 1  || const == 12 || const == 13 || const == 14
Endcost =  - 0.01*mfuel(end) - ThirdStagePayloadMass;
end

EndpointCost = Endcost;



if const == 3 

% RunningCost =((q-50000).^2+1000000)/1000000 + [0 0.005*abs(omegadot)]; % if a cost does not work, try loosening it 
% RunningCost =((q-50000).^2+2000000)/2000000 + [0 0.005*abs(omegadot)];
% RunningCost =((q-50000).^2+2000000)/2000000 + abs(thetadot);
% RunningCost =((q-50000).^2+2000000)/2000000; %looser

% RunningCost =((q-50000).^2+10000)/10000 + [0.005*abs(omegadot)]; % tighter
% RunningCost =((q-50000).^2+50000)/50000 + [0.005*abs(omegadot)]; % tighter, this works as of 1/2/17
RunningCost =((q-50000).^2+100000)/100000 + [0.005*abs(omegadot)]; % tighter
RunningCost =((q-50000).^2+100000)/100000 ; 

% RunningCost =((q-50000).^2+500000)/500000 + [0.005*abs(omegadot)]; % tighter
% RunningCost =((q-50000).^2+5000000)/5000000 + [0.005*abs(omegadot)]; % tighter

% RunningCost =((q-50000).^2+100000)/100000;

elseif const == 1  || const == 12 || const == 13 || const == 14
    
%smoothing functions (can be adjusted depending on needs, remove if not working
% omegadot = diff(thetadot)./diff(time);
%  RunningCost = [0 0.005*abs(omegadot)]; % for smoothing 50kPa
% RunningCost = [0 0.001*abs(omegadot)]; %for smoothing 45kPa and 55kPa and
% high drag

%     RunningCost = Penalty + abs(omegadot);
% RunningCost = Penalty + abs(omegadot) + Fd/1e5 + [0 abs(diff(flapdeflection))]/100;

% RunningCost = Penalty + abs(omegadot) + Fd/1e5;

%  RunningCost = Penalty;
% RunningCost = Penalty + [0 abs(diff(flapdeflection))]/100;
% RunningCost = Penalty + abs(omegadot);
%   RunningCost = Penalty + [0 10*abs(diff(omegadot))];


% IspCost = zeros(1,length(IspNet));
% for i = 1:length(IspNet)
% if IspNet(i) < 0
%     IspCost(i) = IspNet(i)^2;
% end
% end

% qCost = zeros(1,length(q));
% for i = 1:length(q)
% if q(i) < 20000
%     qCost(i) =(q(i)-20000)^2;
% end
% end

% RunningCost =Penalty + 0.1*abs(omegadot) + IspCost + qCost;

%   RunningCost =Penalty*2 + 0.001*abs(omegadot); %including omegadot makes
%   it more volatile..

%     RunningCost = Penalty*2 ; % The Penalty function ensures that it does not go over 50kPa, but still allows it to search that space. 
RunningCost = Penalty*2.3 ;
    %Omegadot cost smooths the trajectory 
RunningCost = 0;
end

global iterative_V
global iterative_t
global iterative_V_f

if rem(iteration,5000) == 0
    q
    
    iterative_V_f(end+1,:) = cumtrapz(time,v.*sin(theta));
    
    
    cla
    hold on
%     Payload_temp
    ThirdStagePayloadMass
%     phi(end)
%     zeta(end)
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
%     plot(iterative_t(end,:),iterative_V_f(end,:),'Color','r')
    
    
    filename = 'testnew51.gif';
    frame = getframe(10);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');

    drawnow
    
    
   
end

end

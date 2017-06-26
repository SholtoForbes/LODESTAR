function [EndpointCost, RunningCost] = SecondStageReturnCost(primal, algorithm)

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

global iteration

iteration = iteration + 1;


V = primal.states(1,:) ; 

v = primal.states(2,:) ; 

 
gamma = primal.states(3,:) ; 

alpha = primal.states(4,:);

zeta = primal.states(5,:);

phi = primal.states(6,:);

xi = primal.states(7,:);

alphadot  = primal.controls(1,:); %

time = primal.nodes;

global Atmosphere
density = interp1(Atmosphere(:,1),Atmosphere(:,4),V);
q = 0.5*density.*v.^2;

% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
end

RunningCost = 0;

% RunningCost = q.^2;

% EndpointCost = time(end);
EndpointCost = 0;
% EndpointCost = H(end);
% EndpointCost = V(end);
EndpointCost = -zeta(end);
global iterative_V
global iterative_t
global iterative_V_f

if rem(iteration,5000) == 0
    
    iterative_V_f(end+1,:) = cumtrapz(time,v.*sin(gamma));
    
    cla
    hold on
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

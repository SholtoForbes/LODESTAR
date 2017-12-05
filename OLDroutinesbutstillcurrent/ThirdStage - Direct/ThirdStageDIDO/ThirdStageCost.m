function [EndpointCost, RunningCost] = ThirdStageCost(primal, algorithm)

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

Alpha = primal.states(4,:) ; 



alphadot  = primal.controls(1,:); %

global time
time = primal.nodes;

phi0 = -0.13;

zeta0 = 1.76;

global rdot
global vdot
global gammadot
global zeta
global phi
global Vec_angle
global inc
global mpayload
global q
global mdot
global m1
global T
global L
global D
global AltF_actual
global Alt1
global Alt_check
global v_check
global gamma_check

[ m1,q,gamma,D,zeta,phi,Vec_angle,T,CL,L, rdot, vdot, gammadot, Alt_check, v_check, gamma_check,mdot] = ThirdStageSim(V, v, gamma, Alpha, phi0, zeta0,nodes,time);


[AltF_actual, vF, Alt, v_forward2, t, mpayload, Alpha, m,gamma_forward2,zeta2,phi2, inc] = ThirdStageSimPostAt(V(end),gamma(end),v(end), phi(end), zeta(end), m1(end)-125.6);



RunningCost = 0;

EndpointCost = -mpayload;

global iterative_V
global iterative_t
global iterative_V_f

if rem(iteration,5000) == 0
    V-Alt_check
    mpayload
    AltF_actual
    rad2deg(Vec_angle)
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

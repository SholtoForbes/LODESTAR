function output = SecondStageEndpoint(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function for Rocket-Scramjet-Rocket System

% This module contains the entire vehicle/system model, as well as defining
% the final cost

% The previous iteration of velocity is initially calculated using primals, and subsequently
% dynamic calculations are performed to produce next velocity step. These
% will eventually converge.

% Global variables are used to pass this velocity to the Dynamics file


%% Import Data %%==========================================================

altF = input.phase.finalstate(1);
vF = input.phase.finalstate(2); 
gammaF = input.phase.finalstate(3); 
% mfuel = input.phase.state(:,4).'; 
% gammadot = input.phase.state(:,5).';
% zeta = input.phase.state(:,6).';

% omegadot  = input.phase.control.'; 
% 
% time = input.phase.time.';

const = input.auxdata.const;

%% THIRD STAGE Payload Matrix %%===========================================
% NEED TO WATCH THIS, IT CAN EXTRAPOLATE BUT IT DOESNT DO IT WELL


if vF > 2850
    ThirdStagePayloadMass = input.auxdata.PayloadGrid(altF,gammaF,vF);
else
    ThirdStagePayloadMass = gaussmf(vF, [300 2850] )*input.auxdata.PayloadGrid(altF,gammaF,2850);
end

% Define Cost =======================================================
% The pseudospectral solver is able to run with no cost at all, this is useful for checking dynamics

if const == 3  || const == 32
Endcost = 0;
% Endcost = 0.001*gamma(end);
% Endcost =1*sum(abs(omegadot));
elseif const == 1  || const == 12 || const == 13 || const == 14  || const == 15
output.objective = -ThirdStagePayloadMass;
end



if const == 3 
    RunningCost =((q-50000).^2+10000)/10000 ; %
elseif const == 32
    RunningCost =((q-50000).^2+10000)/10000;
elseif const == 1  || const == 12 || const == 13 || const == 14  || const == 15
    RunningCost = 0;
end


%% Display Progress %%=====================================================
% iteratio    n = iteration + 1;
% global iterative_V
% global iterative_t
% global iterative_V_f
% if rem(iteration,5000) == 0
%     iterative_V_f(end+1,:) = cumtrapz(time,v.*sin(gamma));
%     cla
%     hold on
%     ThirdStagePayloadMass
%     todisp = ['Total Iterations Done:',num2str(iteration)];
%     disp(todisp);
%     figure(10)
%     iterative_V(end+1,:) = V;
%     iterative_t(end+1,:) = time;
%     if length(iterative_V(:,1)) > 4
%         plot(iterative_t(end-4,:),iterative_V(end-4,:),'Color',[0.8 0.8 0.8])
%     end
%     if length(iterative_V(:,1)) > 2
%         plot(iterative_t(end-2,:),iterative_V(end-2,:),'Color',[0.4 0.4 0.4])
%     end
%     if length(iterative_V(:,1)) > 3
%         plot(iterative_t(end-3,:),iterative_V(end-3,:),'Color',[0.6 0.6 0.6])
%     end
%     plot(iterative_t(end-1,:),iterative_V(end-1,:),'Color',[0.2 0.2 0.2])
%     plot(iterative_t(end,:),iterative_V(end,:),'Color',[0 0 0])
%     filename = 'testnew51.gif';
%     frame = getframe(10);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,filename,'gif','WriteMode','append');
%     drawnow
% end

end

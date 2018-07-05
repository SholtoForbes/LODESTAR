function [eventCost, runningCost] = FirstStageReturnCost(primal)

% global q
Atmosphere = dlmread('atmosphere.txt');


V = primal.states(1,:);	
v = primal.states(2,:);	

alpha = primal.states(4,:);
% density = interp1(Atmosphere(:,1),Atmosphere(:,4),V);
% q = 0.5*density.*v.^2;


eventCost   = 0;
runningCost =0;
% runningCost =abs(alpha);
% runningCost =0.0001*(alpha).^2;
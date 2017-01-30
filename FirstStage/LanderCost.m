function [eventCost, runningCost] = LanderCost(primal)
% global q
Atmosphere = dlmread('atmosphere.txt');


V = primal.states(1,:);	
v = primal.states(2,:);	
alpha = primal.states(5,:);	

density = interp1(Atmosphere(:,1),Atmosphere(:,4),V);
q = 0.5*density.*v.^2;
% eventCost   = -v(end);
% eventCost   = 0;
eventCost   = -v(end)*1000 + (q(end)-55000)^2;
% runningCost =0;
runningCost =abs(alpha);
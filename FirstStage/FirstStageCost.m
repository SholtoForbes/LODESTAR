function [eventCost, runningCost] = LanderCost(primal)

global Atmosphere

V = primal.states(1,:);	
v = primal.states(2,:);	
m = primal.states(3,:);	
gamma = primal.states(4,:);	
alpha = primal.states(5,:);	

density = interp1(Atmosphere(:,1),Atmosphere(:,4),V);
q = 0.5*density.*v.^2;
eventCost   = .01*m(1);

runningCost =0;

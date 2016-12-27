function [eventCost, runningCost] = LanderCost(primal)

v = primal.states(2,:);	
alpha = primal.states(5,:);	

% eventCost   = -v(end);
eventCost   = 0;

% runningCost =0;
runningCost =abs(alpha);
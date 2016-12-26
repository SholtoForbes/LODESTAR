function [eventCost, runningCost] = LanderCost(primal)

v = primal.states(2,:);	
% eventCost   = -v(end);
eventCost   = 0;

runningCost =0;
% That's it!  Remember to fill the first output first!
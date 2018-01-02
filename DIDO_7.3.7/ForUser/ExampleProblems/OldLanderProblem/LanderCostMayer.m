function [eventCost, runningCost] = LanderCost(primal)
% Cost function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eventCost   = -primal.states(3,end);
runningCost = 0;
% That's it!  Remember to fill the first output first!
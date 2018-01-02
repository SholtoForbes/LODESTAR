function [endPointCost, integrandCost] = stageCost(primal)
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Maximize final altitude

endPointCost    = -primal.states(1,end);   
integrandCost   = 0;

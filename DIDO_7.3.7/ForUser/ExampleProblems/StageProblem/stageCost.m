function [endPointCost, integrandCost] = stageCost(primal)
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Example 2, pp.403-404 in I.M.Ross and F. Fahroo,
% "Pseudospectral Knotting Methods for Solving Optimal Control Problems,"
% Journal of Guidance, Control and Dynamics, Vol.27, No.3, May-June 2004.
%-----------------------------------------------
% Maximize final altitude: r = primal.states(1,:)

endPointCost    = -primal.states(1,end);   
integrandCost   = 0;

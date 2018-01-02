function eventConditions = stageEvents(primal)
% Events function for the Staging Problem
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Example 2, pp.403-404 in I.M.Ross and F. Fahroo,
% "Pseudospectral Knotting Methods for Solving Optimal Control Problems,"
% Journal of Guidance, Control and Dynamics, Vol.27, No.3, May-June 2004.
%-----------------------------------------------%
global CONSTANTS

r0 = primal.states(1,1);    rf = primal.states(1,end);
v0 = primal.states(2,1);
m0 = primal.states(3,1);    mf = primal.states(3,end);

%--------------------------------------------------------------------------
% Collect all the left and right limit points
Left    = primal.indices.left;      
Right   = primal.indices.right;
%--------------------------------------------------------------------------
% If this problem had more than two stages, Left and Right would be row
% vectors; here they are just numbers
%------------------------------------------------------------------------
preSeparation_r = primal.states(1, Left);   postSeparation_r = primal.states(1, Right);
preSeparation_v = primal.states(2, Left);   postSeparation_v = primal.states(2, Right);
preSeparation_m = primal.states(3, Left);   postSeparation_m = primal.states(3, Right);
%---------------------------------------------------------------------------------------
% If this problem had more than two stages, preSeparation_r would be a row
% vector of all the left limit points of r.  Ditto for other variables. 
%----------------------------------------------------------------------
%% pre-allocate 8 event conditions
eventConditions = zeros(8,1);
%%
eventConditions(1) = r0;
eventConditions(2) = v0;
eventConditions(3) = m0;

% -- see Equations 63-66 in the above paper.

stage1FuelUsed = m0 - preSeparation_m;
eventConditions(4) = preSeparation_r - postSeparation_r;
eventConditions(5) = preSeparation_v - postSeparation_v;
eventConditions(6) = stage1FuelUsed;
eventConditions(7) = postSeparation_m;

eventConditions(8) = mf;
function eventConditions = stageEvents(primal)
% Events function for the Staging Problem
%--------------------------------------------------------------
% Example file for DIDO 
% TBD DIDO User's Manual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS

r0 = primal.states(1,1);    rf = primal.states(1,end);
v0 = primal.states(2,1);
m0 = primal.states(3,1);    mf = primal.states(3,end);

%-------------------------------------------------------------------------------------
Left    = primal.indices.left;      
Right   = primal.indices.right;

preSeparation_r = primal.states(1, Left);   postSeparation_r = primal.states(1, Right);
preSeparation_v = primal.states(2, Left);   postSeparation_v = primal.states(2, Right);
preSeparation_m = primal.states(3, Left);   postSeparation_m = primal.states(3, Right);
%---------------------------------------------------------------------------------------
%% pre-allocate
eventConditions = zeros(8,1);
%%
eventConditions(1) = r0;
eventConditions(2) = v0;
eventConditions(3) = m0;

stage1FuelUsed = m0 - preSeparation_m;
eventConditions(4) = preSeparation_r - postSeparation_r;
eventConditions(5) = preSeparation_v - postSeparation_v;
eventConditions(6) = stage1FuelUsed;
eventConditions(7) = postSeparation_m;

eventConditions(8) = mf;
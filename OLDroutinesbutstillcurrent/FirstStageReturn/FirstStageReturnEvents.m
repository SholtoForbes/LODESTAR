function eventConditions = FirstStageReturnEvents(primal)
% Events function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO 
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS

h0 = primal.states(1,1);		
v0 = primal.states(2,1);		
gamma0 = primal.states(3,1);		
alpha0 = primal.states(4,1);
zeta0 = primal.states(5,1);

vf = primal.states(2,end);
gammaf = primal.states(3,end);

% eventConditions = zeros(7, 1); % alternatively, write "eventConditions(5) = vf" first.

%===========================================================
eventConditions(1) = h0;
eventConditions(2) = v0;
eventConditions(3) = gamma0;
eventConditions(4) = zeta0;
%-----------------------------------------------------------
eventConditions(5) = vf;
eventConditions(6) = gammaf;
%===========================================================
% all done!

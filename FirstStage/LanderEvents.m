function eventConditions = LanderEvents(primal)
% Events function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO 
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CONSTANTS

h0 = primal.states(1,1);		hf = primal.states(1,end);
v0 = primal.states(2,1);		vf = primal.states(2,end);
m0 = primal.states(3,1);		mf = primal.states(3,end);
gamma0 = primal.states(4,1);		gammaf = primal.states(4,end);
alpha0 = primal.states(5,1);
zeta0 = primal.states(6,1);		zetaf = primal.states(6,end);

% eventConditions = zeros(7, 1); % alternatively, write "eventConditions(5) = vf" first.

%===========================================================
eventConditions(1) = h0;
eventConditions(2) = v0;
eventConditions(3) = m0;
eventConditions(4) = gamma0;
eventConditions(5) = alpha0;
%-----------------------------------------------------------

eventConditions(6) = mf;
eventConditions(7) = zetaf;
% eventConditions(8) = hf;
% eventConditions(8) = gammaf;
%===========================================================
% all done!

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

% Obviously all this reassingment is not necessary but it's easy
% to read the following.  Some penalty in computational performance 
% is incurred.

% preallocate boundary conditions (i.e. event conditions) for good MATLAB computing

eventConditions = zeros(5, 1); % alternatively, write "eventConditions(5) = vf" first.

%===========================================================
eventConditions(1) = h0;
eventConditions(2) = v0;
eventConditions(3) = m0;
%-----------------------------------------------------------
eventConditions(4) = hf;
eventConditions(5) = vf;
%===========================================================
% all done!

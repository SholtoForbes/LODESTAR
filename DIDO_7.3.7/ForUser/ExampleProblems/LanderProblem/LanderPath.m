function ThrustLimits = LanderPath(primal)
% Path constraints function for the Moon-Landing Problem 
%--------------------------------------------------------------
% Example file for DIDO 
% For DIDO User's Manual
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ThrustLimits = primal.primal.states(2,:);

%all done!  Why, were you expecting something else?
function endpointFunction = BadBrac1Events(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint constraint function, e, for the Brac: 1 Problem 
% Template for A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global UNITS

x0Bar = primal.states(1,1);		xfBar = primal.states(1,end);
y0Bar = primal.states(2,1);		yfBar = primal.states(2,end);
v0Bar = primal.states(3,1);		vfBar = primal.states(3,end);



% preallocate the endpointFunction evaluation for good MATLAB computing

endpointFunction = zeros(5,1); % t0 is specified in the main file

%===========================================================
endpointFunction(1) = x0Bar;
endpointFunction(2) = y0Bar;
endpointFunction(3) = v0Bar;

%-----------------------------------------------------------
endpointFunction(4) = xfBar;
endpointFunction(5) = yfBar;  
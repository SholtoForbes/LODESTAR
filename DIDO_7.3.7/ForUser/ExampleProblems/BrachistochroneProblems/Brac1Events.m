function endpointFunction = Brac1Events(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint function for the Brac: 1 Problem 
% Template for a A Beginner's Guide to DIDO 
% I. Michael Ross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = primal.states(1,1);		xf = primal.states(1,end);
y0 = primal.states(2,1);		yf = primal.states(2,end);
v0 = primal.states(3,1);		vf = primal.states(3,end);


% preallocate the endpointFunction evaluation for good MATLAB computing

endpointFunction = zeros(5,1); % t0 is specified in the problem file

%===========================================================
endpointFunction(1) = x0;
endpointFunction(2) = y0;
endpointFunction(3) = v0;

%-----------------------------------------------------------
endpointFunction(4) = xf;
endpointFunction(5) = yf; 
%-----------------------------------------------------------
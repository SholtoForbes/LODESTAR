function endpointFunction = TwoStage2DEvents(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event Function for 2D Problem
% Written by Sholto Forbes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Stage
global const
     

V0 = primal.states(1,1); 
Vf = primal.states(1,end); 
v0 = primal.states(2,1);  
gamma0 = primal.states(3,1);
gammaf = primal.states(3,end); 
mfuel0 = primal.states(4,1);
mfuelf = primal.states(4,end);

omega0 = primal.states(5,1);
omegaf = primal.states(5,end);

zetaf = primal.states(6,end);

if const == 3 || const == 32
endpointFunction = zeros(4,1);
%===========================================================
    endpointFunction(1) = v0;
endpointFunction(2) = mfuel0;
endpointFunction(3) = mfuelf;
endpointFunction(4) = zetaf;
endpointFunction(5) = gammaf;
else
endpointFunction = zeros(5,1);
%===========================================================
    endpointFunction(1) = v0;
endpointFunction(2) = mfuel0;
endpointFunction(3) = mfuelf;
endpointFunction(4) = zetaf;
% endpointFunction(5) = Vf; 
endpointFunction(5) = gammaf;
end
end

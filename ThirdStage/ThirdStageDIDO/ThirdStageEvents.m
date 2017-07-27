function endpointFunction = ThirdStageEvents(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event Function for 2D Problem
% Written by Sholto Forbes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Stage
global const
     global q
     global AltF_actual

V0 = primal.states(1,1); 
Vf = primal.states(1,end); 
v0 = primal.states(2,1);  
gamma0 = primal.states(3,1);  
gammaf = primal.states(3,end);  



endpointFunction = zeros(4,1); % 

endpointFunction(1) = V0;
endpointFunction(2) = v0;
endpointFunction(3) = gamma0;

endpointFunction(4) = AltF_actual;

end

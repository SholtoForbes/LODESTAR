function endpointFunction = SecondStageReturnEvents(primal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event Function for 2D Problem
% Written by Sholto Forbes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Stage
global const
     

V0 = primal.states(1,1); 
Vf = primal.states(1,end); 
v0 = primal.states(2,1);  
vf = primal.states(2,end);  
gamma0 = primal.states(3,1);  
gammaf = primal.states(3,end);
% omega0 = primal.states(4,1);
% omegaf = primal.states(4,end);
zeta0 = primal.states(5,1);
zetaf = primal.states(5,end);
phi0 = primal.states(6,1);
phif = primal.states(6,end);
xi0 = primal.states(7,1);
xif = primal.states(7,end);



endpointFunction = zeros(7,1); % 

endpointFunction(1) = V0;
endpointFunction(2) = v0;
endpointFunction(3) = gamma0;
endpointFunction(4) = zeta0;
endpointFunction(5) = phi0;
endpointFunction(6) = xi0;
% endpointFunction(7) = vf;
% endpointFunction(7) = zetaf;
% endpointFunction(7) = Vf;
endpointFunction(7) = gammaf;
end

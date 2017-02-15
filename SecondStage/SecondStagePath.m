function path = SecondStagePath(primal)
v = primal.states(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global q

DynamicPressure = q;
v0 = v(1);
path = [DynamicPressure ;v0*ones(1,length(v))];





function path = SecondStagePath(primal)
v = primal.states(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global q

DynamicPressure = q;
v0 = v(1);

vfunc = v0 -1524*(1-1e-10*(50000-q(1))^2);
% path = [DynamicPressure ;v0*ones(1,length(v))];
path = [DynamicPressure ;vfunc*ones(1,length(v))];
% path = DynamicPressure ;





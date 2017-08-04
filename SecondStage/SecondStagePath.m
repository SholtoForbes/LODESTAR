function path = SecondStagePath(primal)
V = primal.states(1,:);
v = primal.states(2,:);
theta  = primal.states(3, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global scattered
global q
global const
global M
global Alpha
DynamicPressure = q;
AoA = Alpha;
v0 = v(1);


% vfunc = v0 -scattered.FirstStagev(V(1),theta(1));

% qfunc = DynamicPressure(1) - 50000;
% Mfunc = M(1) - 5;

% vfunc = v0 - griddata(scattered.FirstStageData(:,2),scattered.FirstStageData(:,3),scattered.FirstStageData(:,4),V(1),theta(1),'cubic'); %remember this cant extrapolate
% if isnan(vfunc) == true
%     vfunc = v0 -scattered.FirstStagev(V(1),theta(1));
% end
    
if const == 3
%    path = vfunc*ones(1,length(v)); 
else
% path = [DynamicPressure ;vfunc*ones(1,length(v))];

path = [DynamicPressure; AoA];
% path = [DynamicPressure ;Mfunc*ones(1,length(v));qfunc*ones(1,length(v))];
end
% path = DynamicPressure ;





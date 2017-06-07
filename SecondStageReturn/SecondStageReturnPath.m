function path = SecondStageReturnPath(primal)

global q

DynamicPressure = q;

% vfunc = v0 -scattered.FirstStagev(V(1),theta(1));

% qfunc = DynamicPressure(1) - 50000;
% Mfunc = M(1) - 5;

% vfunc = v0 - griddata(scattered.FirstStageData(:,2),scattered.FirstStageData(:,3),scattered.FirstStageData(:,4),V(1),theta(1),'cubic'); %remember this cant extrapolate
% if isnan(vfunc) == true
%     vfunc = v0 -scattered.FirstStagev(V(1),theta(1));
% end
    

path = [DynamicPressure];

end
% path = DynamicPressure ;





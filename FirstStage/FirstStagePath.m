function qend = LanderPath(primal)
v = primal.states(2,:);
h = primal.states(1,:);

Atmosphere = dlmread('atmosphere.txt');
density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
q = 0.5*density.*v.^2;

qend = q(end);



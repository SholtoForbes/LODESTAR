function Cost = COST(x,y,mFuel,scattered)
phase = 'postpitch';

postpitch0 = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) x(1)];

for i = 1:length(x)
postpitch0(5) = x(i);

tspan = [mFuel/156/length(x)*(i-1) mFuel/156/length(x)*i];

[t_postpitch, postpitch] = ode45(@(t,postpitch) rocketDynamics(postpitch,0,phase,scattered), tspan, postpitch0);

postpitch0 = postpitch(end,:);

end

Cost = postpitch;

end
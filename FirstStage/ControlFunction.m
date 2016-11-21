function Thrust = ControlFunction(t_search,t,u)
Thrust = interp1(t,u,t_search);
end

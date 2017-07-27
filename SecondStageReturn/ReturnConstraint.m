function [constraint,ceq] = ReturnConstraint(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,returncond,num_div)

[cost,phi,t,y,q,xi,y_discrete] = BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,num_div);
% 

% return_constraint = (phi(end) - returncond.phi)^2 + (xi(end) - returncond.xi)^2 - returncond.rad^2;

return_constraint = (phi(end) - returncond.phi)^2 - returncond.rad^2;

q_constraint = max(q) - 50000; % dynamic pressure must be less than 50kPa

crash_constraint = -min(y(:,1)); %Altitude must be greater than 0

end_alt_constraint = y(end,1)-200; % end altitude must be less than 100m

etamin_constraint = max(-(y(:,7) + 1))*1000;

etamax_constraint = max(y(:,7) - 1.5)*1000;

alphamin_constraint = max(-y(:,8))*100;

alphamax_constraint = max(y(:,8) - 8)*100;

% Altitude_0 = controls(num_div*2+2:num_div*3);

% MS_constraint = Altitude_0 - y_discrete(1:end-1);

constraint = [q_constraint crash_constraint etamin_constraint etamax_constraint alphamin_constraint alphamax_constraint];
% constraint = [];
% ceq = [MS_constraint];
ceq = [];
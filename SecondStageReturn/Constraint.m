function [constraint,ceq] = Constraint(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,returncond)

[cost,phi,t,y,q,xi] = BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty);
% 

return_constraint = (phi(end) - returncond.phi)^2 + (xi(end) - returncond.xi)^2 - returncond.rad^2;

q_constraint = max(q) - 50000; % dynamic pressure must be less than 50kPa

alt_constraint = -y(end,1); %Altitude must be greater than 0

constraint = [q_constraint return_constraint alt_constraint];
% constraint = [];
ceq = [];
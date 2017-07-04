function [phi,t,y] = BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty)
%BANKOPT A function for banking optimisation
%   Runs ForwardSimReturn and outputs latitude

% Set angle of attack and roll time histories
alpha_hist = controls(1:length(controls)/2);
eta_hist = controls(length(controls)/2+1:length(controls));

y = [];
t = [];

% Alpha = 7; % aoa (deg)
FlapDeflection = 0;
Iterative_0 = Initial_States; % set the iterative initial states to the global initial states
runtime = 400;
for i = 1:length(eta_hist)
    eta = eta_hist(i); % set roll to the corresponding value in the time series
    alpha = alpha_hist(i); % set roll to the corresponding value in the time series
    [t_iterative, y_iterative] = ode45(@(f_t,f_y) ForwardSimReturn(f_y,alpha,eta,Atmosphere,interp,FlapDeflection,mSPARTAN_empty),[runtime/length(eta_hist)*(i-1) runtime/length(eta_hist)*i],Iterative_0);
    Iterative_0 = y_iterative(end,:);
    
    y = [y;y_iterative];
    t = [t;t_iterative];
end


phi = y(end,2);
end


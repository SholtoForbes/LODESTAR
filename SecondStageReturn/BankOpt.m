function [cost,phi,t,y,q,xi,y_discrete] = BankOpt(controls,Initial_States,Atmosphere,interp,mSPARTAN_empty,num_div)
%BANKOPT A function for banking optimisation
%   Runs ForwardSimReturn and outputs latitude

% Set angle of attack and roll time histories
% alpha_hist = controls(1:num_div);
% eta_hist = controls(num_div+1:2*num_div);

alphadot_hist = controls(1:num_div);
etadot_hist = controls(num_div+1:2*num_div);

% Altitude_0 = controls(num_div*2+2:num_div*3);

y = [];
t = [];
y_discrete = [];

FlapDeflection = 0;
Iterative_0 = Initial_States; % set the iterative initial states to the global initial states
Iterative_0(end+1) = controls(end-1); % add initial roll
Iterative_0(end+1) = controls(end);   % add initial AoA

runtime = controls(num_div*2+1);

for i = 1:num_div
    etadot = etadot_hist(i); % set roll to the corresponding value in the time series
    alphadot = alphadot_hist(i); % set roll to the corresponding value in the time series
    
    [t_iterative, y_iterative] = ode45(@(f_t,f_y) ForwardSimReturn(f_y,alphadot,etadot,Atmosphere,interp,FlapDeflection,mSPARTAN_empty),[runtime/length(etadot_hist)*(i-1) runtime/length(etadot_hist)*i],Iterative_0);
    if i < num_div
        Iterative_0 = y_iterative(end,:);
%         Iterative_0(1) = Altitude_0(i);
    end
    y = [y;y_iterative];
    t = [t;t_iterative];
    
    y_discrete(i) = y(end,1);
end

V = y(:,1);
% min(V)
v = y(:,4);
rho = spline( Atmosphere(:,1),  Atmosphere(:,4), V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

phi = y(:,2);
xi = y(:,6);

cost = 0.01*v(end);
end


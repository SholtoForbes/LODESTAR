%-----------------------------------------------------------------------------
% Ground Distance Calculation, JC May 2017
% Flat earth formula
ground_distance = 6371 * sqrt((phi(end)-phi(1))^2 + (cos(phi(end)-phi(1))*(xi(end)-xi(1)))^2)

% Ellipsoidal Earth formula
K_1 = 111.13209-0.56605*cos(2*(phi(end)-phi(1))) + 0.00120*cos(4*(phi(end)-phi(1)));
K_2 = 111.41513*cos((phi(end)-phi(1))) - 0.09455*cos(3*(phi(end)-phi(1))) + 0.00012*cos(5*(phi(end)-phi(1)));
ground_distance_ellip = sqrt((K_1*radtodeg(phi(end) - phi(1)))^2 + (K_2*radtodeg(xi(end)-xi(1)))^2)

% check by computing ground distances at each point in time history
% forward integrating

ground_distance_list = zeros(1,length(t)-1);

for n = 2:length(t)-1
   ground_distance_list(n) = v(n)*(t(n+1)-t(n))*cos((gamma(n)+gamma(n-1))/2) + ground_distance_list(n-1); 
end

ground_distance_list(end)
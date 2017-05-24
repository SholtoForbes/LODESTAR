function [dfuel, Fueldt, a, q, M, Fd, Thrust, flapdeflection, Alpha, rho,lift,zeta,phi,eq,zetadot] = VehicleModel(time, theta, V, v, mfuel, nodes,scattered, gridded, const,thetadot, Atmosphere, SPARTAN_SCALE,zeta)

% =======================================================
% Vehicle Model
% =======================================================
A = 62.77*SPARTAN_SCALE^(2/3); % reference area (m^2)

eta = .0*ones(1,length(time)); % Roll angle

% eta = 0.3 - 0.0001*time;

%Gravity
g = 9.81;

dt_array = time(2:end)-time(1:end-1); % Time change between each node pt


mstruct = 4910.5 + 3300; % mass of everything but fuel from dawids work

m = mfuel + mstruct;

%===================================================
%
% SECOND STAGE
%
%===================================================

%======================================================

%Rotational Coordinates =================================================
%=================================================

global xi
xi = zeros(1,length(time));
phi = zeros(1,length(time));
% zeta = zeros(1,length(time));

% phi(1) = -0.2138;
phi(1) = -0.264;

% zeta(1) = deg2rad(97);
% zeta(1) = 1.696;
% zeta(1) = 1.699;

r = V + 6371000;
i= 1;


[xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);

for i = 2:length(time)
xi(i) = xi(i-1) + xidot(i-1)*(time(i) - time(i-1));
phi(i) = phi(i-1) + phidot(i-1)*(time(i) - time(i-1));
% zeta(i) = zeta(i-1) + zetadot(i-1)*(time(i) - time(i-1));

[xidot(i),phidot(i),zetadot(i), lift_search(i)] = RotCoords(r(i),xi(i),phi(i),theta(i),v(i),zeta(i),m(i),eta(i), thetadot(i),const);

end


% Aero Data =============================================================
c = spline( Atmosphere(:,1),  Atmosphere(:,5), V); % Calculate speed of sound using atmospheric data

rho = spline( Atmosphere(:,1),  Atmosphere(:,4), V); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

%-heating---------------------------
% kappa = 1.7415e-4;
% Rn = 1; %effective nose radius (m) (need to change this, find actual value)
% 
% heating_rate = kappa*sqrt(rho./Rn).*v.^3; %watts
% 
% Q = zeros(1,length(time));
% Q(1) = 0;
% 
% for i = 1:length(dt_array)
%     Q(i+1) = heating_rate(i)*dt_array(i) + Q(i);
% end

% Control =================================================================

% % determine aerodynamics necessary for trim
% [Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search);
[Fd, Alpha, flapdeflection,lift] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search,c,rho,Atmosphere);



% Alpha
if const == 14
    Fd = 1.1*Fd; % for L/D testing 
end

% THRUST AND MOTION ==================================================================

    
%Fuel Cost ===========================================================================



% kpa50_alt = interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000./v.^2);% altitude at each velocity, if it were a constant 50kPa q trajectory (used to compare with communicator matrix results) 
% kpa50_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), kpa50_alt); % Calculate density using atmospheric data

% Calculate temperature and pressure ratios
% if const == 1 || const == 14
% %     temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), V);
%     
%     Penalty = zeros(1,length(time));
% 
%     for i = 1:length(time)
%         if q(i) < 50000
% 
%         else
% 
%             Penalty(i) = q(i)/50000-1; 
% 
%         end
%     end
% elseif const == 12
% 
% 
%     Penalty = zeros(1,length(time));
% 
%     for i = 1:length(time)
%         if q(i) < 55000
% 
%         else
% 
%             Penalty(i) = q(i)/55000-1; 
% 
%         end
%     end
% elseif const == 13
% 
%     Penalty = zeros(1,length(time));
% 
%     for i = 1:length(time)
%         if q(i) < 45000
% 
%         else
%             Penalty(i) = q(i)/45000-1; 
% 
%         end
%     end
% elseif const == 3 || const == 31
% 
%     Penalty = 0;
% %     t_ratio = temp_actual./kpa50_temp;
% 
% end


lift = lift_search;

T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), V); 
P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), V); 

[Isp,Fueldt,eq] = RESTM12int(M, Alpha, scattered, SPARTAN_SCALE,T0,P0);

for i = 1:length(time)
  if q(i) < 20000
        Isp(i) = Isp(i)*gaussmf(q(i),[1000,20000]);
  end  
end

% for i = 1:length(time)
%   if q(i) < 30000
%         Isp(i) = Isp(i)*gaussmf(q(i),[1000,30000]);
%   end  
% end

Thrust = Isp.*Fueldt*9.81.*cos(deg2rad(Alpha));

fuelchange_array = -Fueldt(1:end-1).*dt_array ;

dfuel = sum(fuelchange_array); %total change in 'fuel' this is negative

v_H = v.*cos(theta);

gravity = m.*(- 6.674e-11.*5.97e24./(V + 6371e3).^2 + v_H.^2./(V + 6371e3)); %Includes Gravity Variation and Centripetal Force 

a = ((Thrust - (Fd + gravity.*sin(theta))) ./ m );

% =========================================================================


end









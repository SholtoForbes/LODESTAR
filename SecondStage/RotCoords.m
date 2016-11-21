function [xidot,phidot,zetadot, lift_search] = RotCoords(r,xi,phi,gamma,v,zeta,m,eta, gammadot,const)
% Determination of motion in rotating coordinates

%xi  Longitude (rad)
%phi  Latitude (rad)
%gamma  Flight Path Angle (rad)
%zeta  Heading Angle (rad)
% 
mu_E = 3.986e14; % m^3/s^2 Earth Gravitational Parameter
omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate

if const == 31
lift_search =(-(v/r - mu_E/(r^2*v))*cos(gamma)*m*v - m*v*cos(phi)*(2*omega_E*cos(zeta) + omega_E^2*r/v*(cos(phi)*cos(gamma)+sin(phi)*sin(gamma)*sin(zeta))))/cos(eta); 
else
lift_search = gammadot*m*v + (-(v/r - mu_E/(r^2*v))*cos(gamma)*m*v - m*v*cos(phi)*(2*omega_E*cos(zeta) + omega_E^2*r/v*(cos(phi)*cos(gamma)+sin(phi)*sin(gamma)*sin(zeta))))/cos(eta); 
end
 
% rdot = v*sin(gamma);

xidot = v*cos(gamma)*cos(zeta)/(r*cos(phi));

phidot = v*cos(gamma)*sin(zeta)/r;

zetadot = lift_search./(m.*v).*sin(eta) + -v/r*tan(phi)*cos(gamma)*cos(zeta) + 2*omega_E*cos(phi)*tan(gamma)*sin(zeta) - omega_E^2*r/(v*cos(gamma))*sin(phi)*cos(phi)*cos(zeta)-2*omega_E*sin(phi);

end


function [rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoords(r,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase)
% Determination of motion in rotating coordinates


mu_E = 3.986e14; % m^3/s^2 Earth Gravitational Parameter
omega_E = 7.292115e-5; % s^-1 Earth Rotation Rate

rdot = v.*sin(gamma);

xidot = v.*cos(gamma).*cos(zeta)./(r.*cos(phi));

phidot = v.*cos(gamma).*sin(zeta)./r;

switch phase
    case 'prepitch'
    gammadot = 0*ones(1,length(gamma)); % Control Trajectory Angle 
    case 'postpitch'
    gammadot = T.*sin(alpha)./(m.*v) + (v./r - mu_E./(r.^2.*v)).*cos(gamma) + L./(m.*v) + cos(phi).*(2.*omega_E.*cos(zeta) + omega_E.^2.*r./v.*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta)));
end


vdot = T.*cos(alpha)./(m) - mu_E.*sin(gamma)./r.^2 -D./m + omega_E.^2.*r.*cos(phi).*(cos(phi).*cos(gamma)+sin(phi).*sin(gamma).*sin(zeta)); 

zetadot = -v./r.*tan(phi).*cos(gamma).*cos(zeta) + 2.*omega_E.*cos(phi).*tan(gamma).*sin(zeta) - omega_E.^2.*r./(v.*cos(gamma)).*sin(phi).*cos(phi).*cos(zeta)-2.*omega_E.*sin(phi);

end


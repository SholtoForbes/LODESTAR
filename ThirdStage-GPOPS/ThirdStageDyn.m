function [rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max] = ThirdStageDyn(r,xi,phi,gamma,v,zeta,m,Alpha,auxdata)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos



time1 = cputime;

r(isnan(r)) = auxdata.Re;

% Atmosphere = dlmread('atmosphere.txt');
% Aero = dlmread('AeroCoeffs.txt');
Atmosphere = auxdata.Atmosphere;

% Drag_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,5));
% 
% Lift_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,6));
% 
% CP_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,7));
% 
% CN_interp = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));
% 
% Max_AoA_interp = scatteredInterpolant(Aero(:,1),Aero(:,4),Aero(:,2));

Drag_interp = auxdata.Drag_interp;

Lift_interp = auxdata.Lift_interp;

CP_interp = auxdata.CP_interp;

CN_interp = auxdata.CN_interp;

Max_AoA_interp = auxdata.Max_AoA_interp;


iteration = 1;


rho_init = spline( Atmosphere(:,1),  Atmosphere(:,4), r(1)-auxdata.Re);
c_init = spline( Atmosphere(:,1),  Atmosphere(:,5), r(1)-auxdata.Re);

q_init = 0.5.*rho_init.*v.^2;
M_init = v./c_init;

%% Set Max AoA
% determine the maximum allowable normal coefficient with a 10 degree limit
% at 50kPa dynamic pressure, and set the max AoA to match this normal
% coefficient
Alt_50 = spline( Atmosphere(:,4),  Atmosphere(:,1), 50000.*2./v(1).^2);
c_50 = spline( Atmosphere(:,1),  Atmosphere(:,5), Alt_50);

M_50 = v(1)./c_50;
% M_50 = 2922.8./c_50;% using a constant velocity

CN_50 = CN_interp(M_50,10);
% CN_50 = CN_interp(M_50,5);
% AoA_max = deg2rad(Max_AoA_interp(M_init,CN_50.*50000./q_init)); %maximum allowable AoA



%Reference area
% A = 0.866; % diameter of 1.05m
A = 0.95; % diameter of 1.1m
g = 9.806; %standard gravity

% the Isp influences the optimal burn mass
% Isp = 437; % from Tom Furgusens Thesis %RL10

Isp = 317.*0.98; %Kestrel, from Falcon 1 users guide, with efficiency reduction

% Isp = 320
% Isp = 446; %HM7B
% Isp = 340; %Aestus 2

%% Define Vehicle Properties
mHS = 130.9; % Heat Shield Mass

% mEng = 100; %RL10
mEng = 52; %Kestrel
% mEng = 165; %HM7B
% mEng = 138; %Aestus 2 ./ RS72 from https:././web.archive.org./web./20141122143945./http:././cs.astrium.eads.net./sp./launcher-propulsion./rocket-engines./aestus-rs72-rocket-engine.html

m(1) = 3300;
% m(1) = 3200;
% m(1) = x(end-2);

% mdot = 14.71; %RL10
% mdot = 9.86977; %Kestrel

mdot = 9.86977.*1.5; %Kestrel Modified

% mdot = 14.72
% mdot = 14.8105; %HM7B
% mdot = 16.5; %Aestus 2



%% Initiate Simulation


% p_spline = spline( Atmosphere(:,1),  Atmosphere(:,3)); % calculate pressure using atmospheric data
% 
% c_spline = spline( Atmosphere(:,1),  Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
% 
% rho_spline = spline( Atmosphere(:,1),  Atmosphere(:,4)); % Calculate density using atmospheric data


atmo_elem = find(r<=85000+auxdata.Re); % Find the elements that are within the atmosphere
exo_elem = find(r>85000+auxdata.Re); % Find the elements that are exoatmospheric


rho(atmo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,4), r(atmo_elem)-auxdata.Re);
c(atmo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,5), r(atmo_elem)-auxdata.Re);
p(atmo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,3), r(atmo_elem)-auxdata.Re);

rho(exo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,4),85000).*gaussmf(r(exo_elem)-auxdata.Re,[100 85000]);
c(exo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,5), 85000);
p(exo_elem) = spline( Atmosphere(:,1),  Atmosphere(:,3),85000).*gaussmf(r(exo_elem)-auxdata.Re,[100 85000]);

rho = rho.';
c = c.';
p = p.';


q(atmo_elem) = 0.5.*rho(atmo_elem).*v(atmo_elem).^2;
M(atmo_elem)  = v(atmo_elem) ./c(atmo_elem) ;
AoA_max(atmo_elem) = deg2rad(Max_AoA_interp(M(atmo_elem),CN_50.*50000./q(atmo_elem))); %maximum allowable AoA

q(exo_elem) = 0.5.*rho(exo_elem).*v(exo_elem).^2;
M(exo_elem)  = v(exo_elem) ./c(exo_elem) ;
AoA_max(exo_elem) = inf; %maximum allowable AoA

q = q.';
M = M.';
AoA_max = AoA_max.';

    
T  = Isp.*mdot.*g - p .*A; % Thrust (N)

 

% AoA_max = 10;
 
    CD  = Drag_interp(M ,rad2deg(Alpha ));
    CL  = Lift_interp(M ,rad2deg(Alpha ));
    CN  = CN_interp(M ,rad2deg(Alpha ));
    cP  = CP_interp(M ,rad2deg(Alpha ));


    D  = 1./2.*rho.*(v.^2).*A.*CD ;
    L  = 1./2.*rho.*(v.^2).*A.*CL ; % Aerodynamic lift
    N  = 1./2.*rho.*(v.^2).*A.*CN ;
    
   
    % Thrust vectoring


        Vec_angle  = asin((cP.*1.1)./2.9554.*N./T ); % calculate the thrust vector angle necessary to resist the lift force moment. cP is in ref lengths



% vec_overmax = find(L  > T*sin(deg2rad(80))) % this is not a limit, it just stops it going imaginary
if any(L  > T*sin(deg2rad(80)))
Vec_angle(L  > T*sin(deg2rad(80)))  = deg2rad(80);
end




% Vec_angle = Vec_angle.'

[rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(r,xi,phi,gamma,v,zeta,L,D,T,m,Alpha,Vec_angle);

end


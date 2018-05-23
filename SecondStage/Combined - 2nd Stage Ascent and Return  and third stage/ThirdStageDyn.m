function [rdot,xidot,phidot,gammadot,vdot,zetadot, mdot, Vec_angle, AoA_max, T, L, D, q] = ThirdStageDyn(alt,gamma,v,m,Alpha,time,auxdata, Alphadot, phi, zeta)
% Function for simulating the Third Stage Rocket Trajectory
% Created by Sholto Forbes-Spyratos



time1 = cputime;

alt(isnan(alt)) = 40000;

Atmosphere = auxdata.interp.Atmosphere;

Drag_interp = auxdata.interp.Drag_interp3;

Lift_interp = auxdata.interp.Lift_interp3;

CP_interp = auxdata.interp.CP_interp3;

CN_interp = auxdata.interp.CN_interp3;

Max_AoA_interp = auxdata.interp.Max_AoA_interp3;


iteration = 1;


rho_init = ppval(auxdata.interp.rho_spline, alt(1));
c_init = ppval(auxdata.interp.c_spline, alt(1));

q_init = 0.5.*rho_init.*v.^2;
M_init = v./c_init;

%% Set Max AoA
% determine the maximum allowable normal coefficient with a 10 degree limit
% at 50kPa dynamic pressure, and set the max AoA to match this normal
% coefficient
Alt_50 = spline( Atmosphere(:,4),  Atmosphere(:,1), 50000.*2./v(1).^2);
c_50 = ppval(auxdata.interp.c_spline, Alt_50);

M_50 = v(1)./c_50;
% M_50 = 2922.8./c_50;% using a constant velocity

CN_50 = CN_interp(M_50,10);
% CN_50 = CN_interp(M_50,5);
% AoA_max = deg2rad(Max_AoA_interp(M_init,CN_50.*50000./q_init)); %maximum allowable AoA



%Reference area

%BASELINE
A = 0.95; % diameter of 1.1m BASELINE

% % Diameter 0.9m
% A = 0.63617;

% Diameter 1m
% A = 0.785;

g = 9.806; %standard gravity

% the Isp influences the optimal burn mass
% Isp = 437; % from Tom Furgusens Thesis %RL10


% BASELINE
Isp = 317.*0.98; %Kestrel, from Falcon 1 users guide, with efficiency reduction. BASELINE

% % Diameter 0.9m
% Isp = 317.*0.96; %Kestrel, from Falcon 1 users guide, with efficiency reduction.

% Diameter 1m
% Isp = 317.*0.97; %Kestrel, from Falcon 1 users guide, with efficiency reduction. 


% Isp = 317;
% Isp = 320
% Isp = 446; %HM7B
% Isp = 340; %Aestus 2

%% Define Vehicle Properties
%BASELINE
mHS = 130.9; % Heat Shield Mass. BASELINE

% %Diameter 0.9m
% mHS = 109.3;

%Diameter 1m
% mHS = 120.341;

% mEng = 100; %RL10
mEng = 52; %Kestrel. BASELINE
% mEng = 165; %HM7B
% mEng = 138; %Aestus 2 ./ RS72 from https:././web.archive.org./web./20141122143945./http:././cs.astrium.eads.net./sp./launcher-propulsion./rocket-engines./aestus-rs72-rocket-engine.html

% m(1) = 3300;
% m(1) = 3200;
% m(1) = x(end-2);

% mdot = 14.71; %RL10
% mdot = 9.86977; %Kestrel

mdot = 9.86977.*1.5; %Kestrel Modified
% mdot = 9.86977*1.2;
% mdot = 14.72
% mdot = 14.8105; %HM7B
% mdot = 16.5; %Aestus 2


% Moment of inertia
% calculated in OneNote, thirdstage, 4 october 2017

%BASELINE
CG = 4.531; % m from end. calculated by hand, but close to the one given from CREO. taken from the end of the rocket. BASELINE, not changed for diameter change (differences are minimal)

% % Diameter 0.9m
% CG =  5.110-0.7521;

%Baseline
% I = 6628; % From Creo

% % Diameter 0.9m
% I = 4864;

% Diameter 1m
I = 5700;

%% Initiate Simulation


atmo_elem = find(alt<=85000); % Find the elements that are within the atmosphere
exo_elem = find(alt>85000); % Find the elements that are exoatmospheric


rho(atmo_elem) = ppval(auxdata.interp.rho_spline, alt(atmo_elem));
c(atmo_elem) = ppval(auxdata.interp.c_spline, alt(atmo_elem));
p(atmo_elem) = ppval(auxdata.interp.p_spline, alt(atmo_elem));

rho(exo_elem) = ppval(auxdata.interp.rho_spline,85000).*gaussmf(alt(exo_elem),[100 85000]);
c(exo_elem) = ppval(auxdata.interp.c_spline, 85000);
p(exo_elem) = ppval(auxdata.interp.p_spline,85000).*gaussmf(alt(exo_elem),[100 85000]);

rho = rho.';
c = c.';
p = p.';


q(atmo_elem) = 0.5.*rho(atmo_elem).*v(atmo_elem).^2;
M(atmo_elem)  = v(atmo_elem) ./c(atmo_elem) ;
AoA_max(atmo_elem) = deg2rad(Max_AoA_interp(M(atmo_elem),CN_50.*50000./q(atmo_elem))); %maximum allowable AoA

q(exo_elem) = 0.5.*rho(exo_elem).*v(exo_elem).^2;
M(exo_elem)  = v(exo_elem) ./c(exo_elem) ;
AoA_max(exo_elem) = deg2rad(30); %maximum allowable AoA

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


        Vec_angle  = asin((cP.*1.1)./(CG-1.5).*(N)./T - Alphadot.*I./((CG-1.5).*T)); % calculate the thrust vector angle necessary to resist the lift force moment. cP is in ref lengths



% vec_overmax = find(L  > T*sin(deg2rad(80))) % this is not a limit, it just stops it going imaginary
if any(L  > T*sin(deg2rad(80)))
Vec_angle(L  > T*sin(deg2rad(80)))  = deg2rad(80);
end


% phi = auxdata.phi0;
% xi = auxdata.xi0;
% zeta = auxdata.zeta0;
% Vec_angle = Vec_angle.'


[rdot,xidot,phidot,gammadot,vdot,zetadot] = RotCoordsRocket(alt+auxdata.Re,0,phi,gamma,v,zeta,L,D,T,m,Alpha,Vec_angle);


end


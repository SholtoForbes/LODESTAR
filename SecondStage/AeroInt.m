function [Cl_spline_EngineOff,Cd_spline_EngineOff,Cl_spline_EngineOn,Cd_spline_EngineOn,flap_spline_EngineOff,flap_spline_EngineOn] = AeroInt(aero,auxdata,T_L,CG_z)
% Sholto 2017
% This function calculates aerodynamic interpolation splines given inputs
% of aerodynamic data matrices.


aero_EngineOff = aero.aero_EngineOff;
aero_EngineOn = aero.aero_EngineOn;
aero_Engine = aero.aero_Engine;
flapaero = aero.flapaero;
Viscousaero_EngineOff = aero.Viscousaero_EngineOff;
Viscousaero_EngineOn = aero.Viscousaero_EngineOn;

%% Aerodynamic Data - Engine off

interp.flap_momentCl_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,3), 'linear', 'nearest');
interp.flap_momentCd_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,4), 'linear', 'nearest');
interp.flap_momentdef_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,2), 'linear', 'nearest');

interp.Cl_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,3));
interp.Cd_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,4));
interp.Cm_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,5));

interp.Cl_scattered_Viscousaero_EngineOff = scatteredInterpolant(Viscousaero_EngineOff(:,1),Viscousaero_EngineOff(:,2),Viscousaero_EngineOff(:,3)/1000,Viscousaero_EngineOff(:,4));
interp.Cd_scattered_Viscousaero_EngineOff = scatteredInterpolant(Viscousaero_EngineOff(:,1),Viscousaero_EngineOff(:,2),Viscousaero_EngineOff(:,3)/1000,Viscousaero_EngineOff(:,5));
    
MList_EngineOff = unique(aero_EngineOff(:,1));
% MList_EngineOff(end+1) = MList_EngineOff(end) + 1; % extrapolate for Mach no slightly

AoAList_engineOff = unique(aero_EngineOff(:,2));

altList_engineOff = unique(Viscousaero_EngineOff(:,3)); % Use engine only case for this

% [Mgrid_EngineOff,AOAgrid_EngineOff] = ndgrid(MList_EngineOff,AoAList_engineOff);

[Mgrid_EngineOff,AOAgrid_EngineOff,altgrid_EngineOff] = ndgrid(MList_EngineOff,AoAList_engineOff,altList_engineOff);

% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

% 

Cl_Grid_EngineOff = [];
Cd_Grid_EngineOff = [];
Cm_Grid_EngineOff = [];
flap_Grid = [];


for i = 1:numel(Mgrid_EngineOff)
    M_temp = Mgrid_EngineOff(i);
    AoA_temp = AOAgrid_EngineOff(i);
    alt_temp = altgrid_EngineOff(i);
    
    Cl_temp_EngineOff = interp.Cl_scattered_EngineOff(M_temp,AoA_temp); % Determine coefficients without flap deflections. 
    Cd_temp_EngineOff = interp.Cd_scattered_EngineOff(M_temp,AoA_temp);
    Cm_temp_EngineOff = interp.Cm_scattered_EngineOff(M_temp,AoA_temp);
    
    Cl_temp_ViscousEngineOff = interp.Cl_scattered_Viscousaero_EngineOff(M_temp,AoA_temp,alt_temp);
    Cd_temp_ViscousEngineOff = interp.Cd_scattered_Viscousaero_EngineOff(M_temp,AoA_temp,alt_temp);
    
    
    %determine Flap Component
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOff(M_temp,0); % Determine coefficients with no flap deflection as reference.
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOff(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOff(M_temp,0);
    
    Cl_AoA0_withflaps_temp_EngineOff = interp.flap_momentCl_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0));  % Interpolate for coefficients with flap deflections. Removing portion of SPARTAN moment with no flap deflection.
    Cd_AoA0_withflaps_temp_EngineOff = interp.flap_momentCd_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0));
    
    
    flap_Cl_temp_EngineOff = Cl_AoA0_withflaps_temp_EngineOff - Cl_temp_AoA0; % Remove portion of coefficient caused by SPARTAN. 
    flap_Cd_temp_EngineOff = Cd_AoA0_withflaps_temp_EngineOff - Cd_temp_AoA0;
    %
    
  % Create Grids
    I = cell(1, ndims(Mgrid_EngineOff)); 
    [I{:}] = ind2sub(size(Mgrid_EngineOff),i);
    
    Cl_Grid_EngineOff(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOff+flap_Cl_temp_EngineOff+Cl_temp_ViscousEngineOff; % Add flap deflection coefficients in. This assumes that flap deflection will have equal effect over the range of Mach no.s.
    Cd_Grid_EngineOff(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOff+flap_Cd_temp_EngineOff+Cd_temp_ViscousEngineOff;
    Cm_Grid_EngineOff(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOff;

    flap_Grid_EngineOff(I{(1)},I{(2)},I{(3)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0)) ;
    
    Cl_Grid_test_EngineOff(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOff; % Baseline test case without flap deflection. 
    Cd_Grid_test_EngineOff(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOff;
    Cm_Grid_test_EngineOff(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOff;
    
%     Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp;
end
Cl_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,altgrid_EngineOff,Cl_Grid_EngineOff,'spline','linear');
Cd_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,altgrid_EngineOff,Cd_Grid_EngineOff,'spline','linear');
flap_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,altgrid_EngineOff,flap_Grid_EngineOff,'spline','linear');
Cm_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,altgrid_EngineOff,Cm_Grid_EngineOff,'spline','linear');


%% Aerodynamic Data - Engine on 


interp.Cl_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3));
interp.Cd_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,4));
interp.Cm_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,5));

interp.Cl_scattered_Engine = scatteredInterpolant(aero_Engine(:,1),aero_Engine(:,2),aero_Engine(:,3),aero_Engine(:,4));
interp.Cd_scattered_Engine = scatteredInterpolant(aero_Engine(:,1),aero_Engine(:,2),aero_Engine(:,3),aero_Engine(:,5));
interp.Cm_scattered_Engine = scatteredInterpolant(aero_Engine(:,1),aero_Engine(:,2),aero_Engine(:,3),aero_Engine(:,6));

interp.Cl_scattered_Viscousaero_EngineOn = scatteredInterpolant(Viscousaero_EngineOn(:,1),Viscousaero_EngineOn(:,2),Viscousaero_EngineOn(:,3)/1000,Viscousaero_EngineOn(:,4));
interp.Cd_scattered_Viscousaero_EngineOn = scatteredInterpolant(Viscousaero_EngineOn(:,1),Viscousaero_EngineOn(:,2),Viscousaero_EngineOn(:,3)/1000,Viscousaero_EngineOn(:,5));

MList_EngineOn = unique(aero_EngineOn(:,1));
% MList_EngineOn(end+1) = MList_EngineOn(end) + 1; % extrapolate for Mach no slightly

AoAList_engineOn = unique(aero_EngineOn(:,2));

altList_engineOn = unique(aero_Engine(:,3)); % Use engine only case for this

[Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn] = ndgrid(MList_EngineOn,AoAList_engineOn,altList_engineOn);


% [Mgrid_EngineOn,AOAgrid_EngineOn] = ndgrid(MList_EngineOn,AoAList_engineOn);

Cl_Grid_EngineOn = [];
Cd_Grid_EngineOn = [];
Cm_Grid_EngineOn = [];
flap_Grid = [];

for i = 1:numel(Mgrid_EngineOn)
    M_temp = Mgrid_EngineOn(i);
    AoA_temp = AOAgrid_EngineOn(i);
    alt_temp = altgrid_EngineOn(i);
    
    
    %% PUT THIS BACK IN WHEN ADDING ENGINE, just taken out because of lack of alt right now 29/6
    % Also modify flaps below, and change interp in cehicle model
    
    %Calculate Thrust
%     
    L_ref = 2294.0;
    
    c = ppval(auxdata.interp.c_spline,alt_temp); % Calculate speed of sound using atmospheric data

    v = M_temp*c;

    rho = ppval(auxdata.interp.rho_spline,alt_temp); % Calculate density using atmospheric data

    q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

    T0 = ppval(auxdata.interp.T0_spline, alt_temp); 

    P0 = ppval(auxdata.interp.P0_spline, alt_temp);

    [Isp,Fueldt,eq,q1] = RESTint(M_temp, AoA_temp, auxdata,T0,P0);

    T = Isp.*Fueldt*9.81; % Thrust in direction of motion

    Tm = T*(CG_z-T_L)/(q*L_ref); %Thrust moment coefficient in same reference at the moments generated by CART3D


    %%
    
    
    Cl_temp_Engine = interp.Cl_scattered_Engine(M_temp,AoA_temp,alt_temp);
    Cd_temp_Engine = interp.Cd_scattered_Engine(M_temp,AoA_temp,alt_temp);
    Cm_temp_Engine = interp.Cm_scattered_Engine(M_temp,AoA_temp,alt_temp);

    Cl_temp_EngineOn = interp.Cl_scattered_EngineOn(M_temp,AoA_temp);
    Cd_temp_EngineOn = interp.Cd_scattered_EngineOn(M_temp,AoA_temp);
    Cm_temp_EngineOn = interp.Cm_scattered_EngineOn(M_temp,AoA_temp);
    
    Cl_temp_ViscousEngineOn = interp.Cl_scattered_Viscousaero_EngineOn(M_temp,AoA_temp,alt_temp);
    Cd_temp_ViscousEngineOn = interp.Cd_scattered_Viscousaero_EngineOn(M_temp,AoA_temp,alt_temp);
    
    %determine Flap Component
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOff(M_temp,0); % engine off case used as reference for flaps (which also have engine off)
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOff(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOff(M_temp,0);
    
    Cl_AoA0_withflaps_temp_EngineOn = interp.flap_momentCl_scattered(M_temp,-(Cm_temp_EngineOn+Cm_temp_Engine+Tm-Cm_temp_AoA0));
    Cd_AoA0_withflaps_temp_EngineOn = interp.flap_momentCd_scattered(M_temp,-(Cm_temp_EngineOn+Cm_temp_Engine+Tm-Cm_temp_AoA0)) ;
%     Cl_AoA0_withflaps_temp_EngineOn = interp.flap_momentCl_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0));
%     Cd_AoA0_withflaps_temp_EngineOn = interp.flap_momentCd_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0)) ;
    
    
    flap_Cl_temp_EngineOn = Cl_AoA0_withflaps_temp_EngineOn - Cl_temp_AoA0;
    flap_Cd_temp_EngineOn = Cd_AoA0_withflaps_temp_EngineOn - Cd_temp_AoA0;
    %
    
    % Create Grids
    I = cell(1, ndims(Mgrid_EngineOn)); 
    [I{:}] = ind2sub(size(Mgrid_EngineOn),i);
    
    Cl_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOn+Cl_temp_Engine+Cl_temp_ViscousEngineOn+flap_Cl_temp_EngineOn;
    Cd_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOn+Cd_temp_Engine+Cd_temp_ViscousEngineOn+flap_Cd_temp_EngineOn;
    Cm_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOn;

%    Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp_EngineOn+flap_Cl_temp_EngineOn;
%     Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp_EngineOn+flap_Cd_temp_EngineOn;
%     Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp_EngineOn;

    flap_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOn+Tm-Cm_temp_AoA0)) ;
% flap_Grid_EngineOn(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0)) ;    


%     Cl_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOn;
%     Cd_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOn;
%     Cm_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOn;
%     
%     Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp;
end
Cl_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cl_Grid_EngineOn,'spline','linear');
Cd_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cd_Grid_EngineOn,'spline','linear');
flap_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,flap_Grid_EngineOn,'spline','linear');
Cm_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cm_Grid_EngineOn,'spline','linear');

% Cl_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,Cl_Grid_EngineOn,'spline','linear');
% Cd_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,Cd_Grid_EngineOn,'spline','linear');
% flap_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,flap_Grid_EngineOn,'spline','linear');
% Cm_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,Cm_Grid_EngineOn,'spline','linear');

plotaero = 'no';
if strcmp(plotaero,'yes')
figure(401)
contourf(Mgrid_EngineOff,AOAgrid_EngineOff,Cl_Grid_test_EngineOff,5000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Lift COefficient, No Flaps, Engine Off')
c = colorbar;
c.Label.String = 'Lift Coefficient';


figure(402)
contourf(Mgrid_EngineOff,AOAgrid_EngineOff,Cd_Grid_test_EngineOff,5000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Drag COefficient, No Flaps, Engine Off')
c = colorbar;
c.Label.String = 'Drag Coefficient';


figure(403)
contourf(Mgrid_EngineOff,AOAgrid_EngineOff,Cl_Grid_test_EngineOff./Cd_Grid_test_EngineOff,5000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('L/D, No Flaps, Engine Off')
c = colorbar;
c.Label.String = 'L/D';
shading interp

figure(404)
contourf(Mgrid_EngineOn,AOAgrid_EngineOn,Cl_Grid_test_EngineOn,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Lift COefficient, No Flaps, Engine On')
c = colorbar;
c.Label.String = 'Lift Coefficient';


figure(405)
contourf(Mgrid_EngineOn,AOAgrid_EngineOn,Cd_Grid_test_EngineOn,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Drag Coefficient, No Flaps, Engine On')
c = colorbar;
c.Label.String = 'Drag Coefficient';


figure(406)
contourf(Mgrid_EngineOn,AOAgrid_EngineOn,Cl_Grid_test_EngineOn./Cd_Grid_test_EngineOn,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('L/D, No Flaps, Engine On')
c = colorbar;
c.Label.String = 'L/D';
shading interp

figure(407)
contourf(Mgrid_EngineOff,AOAgrid_EngineOff,Cl_Grid_EngineOff./Cd_Grid_EngineOff,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('L/D, Flaps, Engine Off')
c = colorbar;
c.Label.String = 'L/D';
shading interp
end

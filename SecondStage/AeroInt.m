function [Cl_spline_EngineOff,Cd_spline_EngineOff,Cl_spline_EngineOn,Cd_spline_EngineOn,flap_spline_EngineOff,flap_spline_EngineOn] = AeroInt(aero_EngineOff,aero_EngineOn,flapaero)
% Sholto 2017
% This function calculates aerodynamic interpolation splines given inputs
% of aerodynamic data matrices.

% addpath ./CG14.5
% 
% aero_EngineOff = importdata('SPARTANaero14.5');
% 
% flapaero = importdata('SPARTANaeroFlaps14.5');
% 
% aero_EngineOn = importdata('SPARTANaeroEngineOn14.5');


%% Aerodynamic Data - Engine off

interp.flap_momentCl_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,3), 'linear', 'nearest');
interp.flap_momentCd_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,4), 'linear', 'nearest');
interp.flap_momentdef_scattered = scatteredInterpolant(flapaero(:,1),flapaero(:,5),flapaero(:,2), 'linear', 'nearest');

interp.Cl_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,3));
interp.Cd_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,4));
interp.Cm_scattered_EngineOff = scatteredInterpolant(aero_EngineOff(:,1),aero_EngineOff(:,2),aero_EngineOff(:,5));


MList_EngineOff = unique(aero_EngineOff(:,1));
% MList_EngineOff(end+1) = MList_EngineOff(end) + 1; % extrapolate for Mach no slightly

AoAList_engineOff = unique(aero_EngineOff(:,2));

[Mgrid_EngineOff,AOAgrid_EngineOff] = ndgrid(MList_EngineOff,AoAList_engineOff);

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
    
    Cl_temp_EngineOff = interp.Cl_scattered_EngineOff(M_temp,AoA_temp); % Determine coefficients without flap deflections. 
    Cd_temp_EngineOff = interp.Cd_scattered_EngineOff(M_temp,AoA_temp);
    Cm_temp_EngineOff = interp.Cm_scattered_EngineOff(M_temp,AoA_temp);
    
    
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
    
    Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp_EngineOff+flap_Cl_temp_EngineOff; % Add flap deflection coefficients in. This assumes that flap deflection will have equal effect over the range of Mach no.s.
    Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp_EngineOff+flap_Cd_temp_EngineOff;
    Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp_EngineOff;

    flap_Grid_EngineOff(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0)) ;
    
    Cl_Grid_test_EngineOff(I{(1)},I{(2)}) = Cl_temp_EngineOff; % Baseline test case without flap deflection. 
    Cd_Grid_test_EngineOff(I{(1)},I{(2)}) = Cd_temp_EngineOff;
    Cm_Grid_test_EngineOff(I{(1)},I{(2)}) = Cm_temp_EngineOff;
    
%     Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp;
end
Cl_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,Cl_Grid_EngineOff,'spline','linear');
Cd_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,Cd_Grid_EngineOff,'spline','linear');
flap_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,flap_Grid_EngineOff,'spline','linear');
Cm_spline_EngineOff = griddedInterpolant(Mgrid_EngineOff,AOAgrid_EngineOff,Cm_Grid_EngineOff,'spline','linear');


%% Aerodynamic Data - Engine on 


interp.Cl_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3),aero_EngineOn(:,4));
interp.Cd_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3),aero_EngineOn(:,5));
interp.Cm_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3),aero_EngineOn(:,6));

MList_EngineOn = unique(aero_EngineOn(:,1));
% MList_EngineOn(end+1) = MList_EngineOn(end) + 1; % extrapolate for Mach no slightly

AoAList_engineOn = unique(aero_EngineOn(:,2));

altList_engineOn = unique(aero_EngineOn(:,3));

[Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn] = ndgrid(MList_EngineOn,AoAList_engineOn,altList_engineOn);


Cl_Grid_EngineOn = [];
Cd_Grid_EngineOn = [];
Cm_Grid_EngineOn = [];
flap_Grid = [];

for i = 1:numel(Mgrid_EngineOn)
    M_temp = Mgrid_EngineOn(i);
    AoA_temp = AOAgrid_EngineOn(i);
    alt_temp = altgrid_EngineOn(i);
    
    Cl_temp_EngineOn = interp.Cl_scattered_EngineOn(M_temp,AoA_temp,alt_temp);
    Cd_temp_EngineOn = interp.Cd_scattered_EngineOn(M_temp,AoA_temp,alt_temp);
    Cm_temp_EngineOn = interp.Cm_scattered_EngineOn(M_temp,AoA_temp,alt_temp);
    
    %determine Flap Component
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOff(M_temp,0); % engine off case used as reference for flaps (which also have engine off)
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOff(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOff(M_temp,0);
    
    Cl_AoA0_withflaps_temp_EngineOn = interp.flap_momentCl_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0));
    Cd_AoA0_withflaps_temp_EngineOn = interp.flap_momentCd_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0)) ;
    
    flap_Cl_temp_EngineOn = Cl_AoA0_withflaps_temp_EngineOn - Cl_temp_AoA0;
    flap_Cd_temp_EngineOn = Cd_AoA0_withflaps_temp_EngineOn - Cd_temp_AoA0;
    %
    
    % Create Grids
    I = cell(1, ndims(Mgrid_EngineOn)); 
    [I{:}] = ind2sub(size(Mgrid_EngineOn),i);
    
    Cl_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOn+flap_Cl_temp_EngineOn;
    Cd_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOn+flap_Cd_temp_EngineOn;
    Cm_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOn;

    flap_Grid_EngineOn(I{(1)},I{(2)},I{(3)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0)) ;
    
    Cl_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cl_temp_EngineOn;
    Cd_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cd_temp_EngineOn;
    Cm_Grid_test_EngineOn(I{(1)},I{(2)},I{(3)}) = Cm_temp_EngineOn;
    
%     Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp;
end
Cl_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cl_Grid_EngineOn,'spline','linear');
Cd_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cd_Grid_EngineOn,'spline','linear');
flap_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,flap_Grid_EngineOn,'spline','linear');
Cm_spline_EngineOn = griddedInterpolant(Mgrid_EngineOn,AOAgrid_EngineOn,altgrid_EngineOn,Cm_Grid_EngineOn,'spline','linear');



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

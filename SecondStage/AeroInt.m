function auxdata = AeroInt(aero_EngineOff,aero_EngineOn,flapaero)

% Calculates aerodynamic interpolation splines given inputs of aerodynamic
% data

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

[MList_EngineOff,AOAList_EngineOff] = ndgrid(unique(aero_EngineOff(:,1)),unique(aero_EngineOff(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid_EngineOff = [];
Cd_Grid_EngineOff = [];
Cm_Grid_EngineOff = [];
flap_Grid = [];


for i = 1:numel(MList_EngineOff)
    M_temp = MList_EngineOff(i);
    AoA_temp = AOAList_EngineOff(i);
    
    Cl_temp_EngineOff = interp.Cl_scattered_EngineOff(M_temp,AoA_temp);
    Cd_temp_EngineOff = interp.Cd_scattered_EngineOff(M_temp,AoA_temp);
    Cm_temp_EngineOff = interp.Cm_scattered_EngineOff(M_temp,AoA_temp);
    
    
    %determine Flap Component
    Cd_temp_AoA0 = interp.Cd_scattered_EngineOff(M_temp,0);
    Cl_temp_AoA0 = interp.Cl_scattered_EngineOff(M_temp,0);
    Cm_temp_AoA0 = interp.Cm_scattered_EngineOff(M_temp,0);
    
    Cl_AoA0_withflaps_temp_EngineOff = interp.flap_momentCl_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0));
    Cd_AoA0_withflaps_temp_EngineOff = interp.flap_momentCd_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0)) ;
    
    
    flap_Cl_temp_EngineOff = Cl_AoA0_withflaps_temp_EngineOff - Cl_temp_AoA0;
    flap_Cd_temp_EngineOff = Cd_AoA0_withflaps_temp_EngineOff - Cd_temp_AoA0;
    %
    
  % Create Grids
    I = cell(1, ndims(MList_EngineOff)); 
    [I{:}] = ind2sub(size(MList_EngineOff),i);
    
    Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp_EngineOff+flap_Cl_temp_EngineOff;
    Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp_EngineOff+flap_Cd_temp_EngineOff;
    Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp_EngineOff;

    flap_Grid_EngineOff(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOff-Cm_temp_AoA0)) ;
    
    Cl_Grid_test_EngineOff(I{(1)},I{(2)}) = Cl_temp_EngineOff;
    Cd_Grid_test_EngineOff(I{(1)},I{(2)}) = Cd_temp_EngineOff;
    Cm_Grid_test_EngineOff(I{(1)},I{(2)}) = Cm_temp_EngineOff;
    
%     Cl_Grid_EngineOff(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOff(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOff(I{(1)},I{(2)}) = Cm_temp;
end
auxdata.interp.Cl_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cl_Grid_EngineOff,'spline','linear');
auxdata.interp.Cd_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cd_Grid_EngineOff,'spline','linear');
auxdata.interp.Cm_spline_EngineOff = griddedInterpolant(MList_EngineOff,AOAList_EngineOff,Cm_Grid_EngineOff,'spline','nearest');

%% Aerodynamic Data - Engine on 


interp.Cl_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,3));
interp.Cd_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,4));
interp.Cm_scattered_EngineOn = scatteredInterpolant(aero_EngineOn(:,1),aero_EngineOn(:,2),aero_EngineOn(:,5));

[MList_EngineOn,AOAList_EngineOn] = ndgrid(unique(aero_EngineOn(:,1)),unique(aero_EngineOn(:,2)));
% Cl_Grid = reshape(aero(:,3),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';
% Cd_Grid = reshape(aero(:,4),[length(unique(aero(:,2))),length(unique(aero(:,1)))]).';

Cl_Grid_EngineOn = [];
Cd_Grid_EngineOn = [];
Cm_Grid_EngineOn = [];
flap_Grid = [];

for i = 1:numel(MList_EngineOn)
    M_temp = MList_EngineOn(i);
    AoA_temp = AOAList_EngineOn(i);
    
    Cl_temp_EngineOn = interp.Cl_scattered_EngineOn(M_temp,AoA_temp);
    Cd_temp_EngineOn = interp.Cd_scattered_EngineOn(M_temp,AoA_temp);
    Cm_temp_EngineOn = interp.Cm_scattered_EngineOn(M_temp,AoA_temp);
    
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
    I = cell(1, ndims(MList_EngineOn)); 
    [I{:}] = ind2sub(size(MList_EngineOn),i);
    
    Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp_EngineOn+flap_Cl_temp_EngineOn;
    Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp_EngineOn+flap_Cd_temp_EngineOn;
    Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp_EngineOn;

    flap_Grid_EngineOn(I{(1)},I{(2)}) = interp.flap_momentdef_scattered(M_temp,-(Cm_temp_EngineOn-Cm_temp_AoA0)) ;
    
    Cl_Grid_test_EngineOn(I{(1)},I{(2)}) = Cl_temp_EngineOn;
    Cd_Grid_test_EngineOn(I{(1)},I{(2)}) = Cd_temp_EngineOn;
    Cm_Grid_test_EngineOn(I{(1)},I{(2)}) = Cm_temp_EngineOn;
    
%     Cl_Grid_EngineOn(I{(1)},I{(2)}) = Cl_temp;
%     Cd_Grid_EngineOn(I{(1)},I{(2)}) = Cd_temp;
%     Cm_Grid_EngineOn(I{(1)},I{(2)}) = Cm_temp;
end
auxdata.interp.Cl_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cl_Grid_EngineOn,'spline','linear');
auxdata.interp.Cd_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cd_Grid_EngineOn,'spline','linear');
auxdata.interp.Cm_spline_EngineOn = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Cm_Grid_EngineOn,'spline','nearest');




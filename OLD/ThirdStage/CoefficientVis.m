clear all


aero = dlmread('AeroCoeffs.txt');

Cl_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,6));
Cd_scattered = scatteredInterpolant(aero(:,1),aero(:,2),aero(:,5));

MList  = unique(aero(:,1));
% MList (end+1) = MList (end) + 1; % extrapolate for Mach no slightly

AoAList  = unique(aero(:,2));


[Mgrid ,AOAgrid ] = ndgrid(MList ,AoAList );


Cl_Grid = [];
Cd_Grid_ = [];



for i = 1:numel(Mgrid )
    M_temp = Mgrid (i);
    AoA_temp = AOAgrid (i);
    
    Cl_temp = Cl_scattered (M_temp,AoA_temp); % Determine coefficients without flap deflections. 
    Cd_temp  = Cd_scattered (M_temp,AoA_temp);
   
        I = cell(1, ndims(Mgrid)); 
    [I{:}] = ind2sub(size(Mgrid),i);
   
    Cl_Grid(I{(1)},I{(2)}) = Cl_temp ; % Baseline test case without flap deflection. 
    Cd_Grid(I{(1)},I{(2)}) = Cd_temp ;


end

figure(404)
contourf(Mgrid,AOAgrid,Cl_Grid,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Lift Coefficient')
c = colorbar;
c.Label.String = 'Lift Coefficient';


figure(405)
contourf(Mgrid,AOAgrid,Cd_Grid,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('Drag Coefficient')
c = colorbar;
c.Label.String = 'Drag Coefficient';


figure(406)
contourf(Mgrid,AOAgrid,Cl_Grid./Cd_Grid,3000,'LineWidth',0.)
xlabel('Mach no.')
ylabel('Angle of Attack (deg)')
title('L/D')
c = colorbar;
c.Label.String = 'L/D';
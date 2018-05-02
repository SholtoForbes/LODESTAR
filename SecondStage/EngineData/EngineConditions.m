% Routine to calculate nozzle exit conditions 
% Uses conical shock calculator, the  interpolates
clear all

M0 = 8;
Alpha = 4;
T0 = 220.650
P0 = 2930.49
rho0 = 0.0462674




auxdata.interp.engine_data = dlmread('ENGINEDATA.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.interp.engine_data;


shockdata = dlmread('ShockMat');
[MList_EngineOff,AOAList_EngineOn] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
auxdata.interp.M1gridded = griddedInterpolant(MList_EngineOff,AOAList_EngineOn,M1_Grid,'spline','linear');
auxdata.interp.presgridded = griddedInterpolant(MList_EngineOff,AOAList_EngineOn,pres_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList_EngineOff,AOAList_EngineOn,temp_Grid,'spline','linear');



T1 = auxdata.interp.tempgridded(M0,Alpha).*T0;
P1 = auxdata.interp.presgridded(M0,Alpha).*P0; % note this is at 50kPa, modified by efficiency
M1 = auxdata.interp.M1gridded(M0, Alpha);



Me_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,5));
Pe_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,6));
Te_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,7));
Ge_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,8));
re_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,9));
P1nom_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,10));


Me = Me_interpolator(M1,T1)
Pe = Pe_interpolator(M1,T1)*P1/P1nom_interpolator(M1,T1)
Te = Te_interpolator(M1,T1)
Ge = Ge_interpolator(M1,T1);
re = re_interpolator(M1,T1);





R = 8.314;


rhoe = Pe/re/Te

A = 1.44704*4 % Calculated using creo
mdot = A*Pe/sqrt(Te) * sqrt(Ge/R) * Me * (1+(Ge-1)/2*Me^2)^-((Ge+1)/(2*(Ge-1)))



%% Nondimensionalisation for CART3D

Pe/P0



Te/T0


%%
%For SurfBC

Pe_ = Pe/Ge/P0

rhoe_ = rhoe/rho0

Ue_ = Me*sqrt(Ge*Pe_/rhoe_)

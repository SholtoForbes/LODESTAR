% Routine to calculate nozzle exit conditions 
% Uses conical shock calculator, the  interpolates
clear all

M0 = 5;
Alpha = 0;
alt = 30;

Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 


T0 = ppval(T0_spline,alt*1000)
P0 = ppval(P0_spline,alt*1000)
rho0 = ppval(rho_spline,alt*1000)





% 24km
% T0 = 220.650
% P0 = 2930.49
% rho0 = 0.0462674

%30km
% T0 = 226.650
% P0 = 1171.87
% rho0 = 0.0180119

% 36km
% T0 = 239.850
% P0 = 484.317
% rho0 = 0.00703441

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



T1 = auxdata.interp.tempgridded(M0,Alpha).*T0
P1 = auxdata.interp.presgridded(M0,Alpha).*P0; % note this is at 50kPa, modified by efficiency
M1 = auxdata.interp.M1gridded(M0, Alpha);

% M1 = 5.58210
% T1 = 252.834

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



load gridIsp_eng
grid.Isp_eng = gridIsp_eng;

% gridIsp_eng may have sections at which the Isp is 0. The following finds
% these, and fills them in with linearly intepolated values.
Isp_interpolator = scatteredInterpolant(engine_data(:,1),engine_data(:,2),engine_data(:,3));

for i = 1:30 % must match engineint.m
    for j= 1:30
        % grid.Isp_eng(i,j) = polyvaln(p,[grid.Mgrid_eng(i,j) grid.T_eng(i,j)]);
        if any(grid.Isp_eng(i,j)) == false
            grid.Isp_eng(i,j) = Isp_interpolator(grid.Mgrid_eng(i,j), grid.T_eng(i,j)); % Linearly extrapolate for any data point which the Bivar spline could not solve
        end
    end
end

auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');

%% Equivalence Ratio %%==========================================================
% Import engine data
auxdata.interp.engine_data = dlmread('ENGINEDATA.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.interp.engine_data;

% Create uniform grid of Mach no. and temperature values. 
M_englist = unique(sort(engine_data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(engine_data(:,1)));

T_englist = unique(sort(engine_data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = unique(sort(engine_data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

% Set the equivalence ratio interpolation region %-------------------------
% VERY IMPORTANT

% The interpolators have trouble with equivalence ratio because its equal
% to 1 over a certain Mach no. (causes error in interpolator, as the
% interpolator will find values of equivalence ratio < 1 where they should
% not exist)

% This makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)

% the the maximum of this to around where equivalence ratio stops changing,
% and check the end results

eq_data = [];
j=1;
for i = 1: length(engine_data(:,1))
    if engine_data(i,1) < 5.
        eq_data(j,:) = engine_data(i,:);
        j=j+1;
    end
end

auxdata.interp.equivalence = scatteredInterpolant(eq_data(:,1),eq_data(:,2),eq_data(:,4), 'linear');
grid.eq_eng = auxdata.interp.equivalence(grid.Mgrid_eng,grid.T_eng);
auxdata.interp.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');


[Isp,wf,eq,q1,w] = RESTint(M0, Alpha, auxdata,T0,P0);
Isp
w+wf

T = Isp.*wf*9.81;


R = 8.314;


rhoe = Pe/re/Te;

A = 0.5683*4;

mm = 0.029; % molar mass air

% Tt = Te*(1+(Ge-1)/2*Me^2)
% Pt = Pe + 0.5*rhoe*(Me*sqrt(Ge*re*Te))^2
% 
% mdot = A*Pt/sqrt(Tt) * sqrt(Ge/(re*mm)) * Me * (1+(Ge-1)/2*Me^2)^-((Ge+1)/(2*(Ge-1)))

ve = Me*sqrt(Ge*re*Te);

mdot = rhoe*A*ve

%% Nondimensionalisation for CART3D
%%
%For SurfBC

Pe_ = Pe/1.4/P0

rhoe_ = rhoe/rho0

Ue_ = sqrt(Ge/1.4*(Me*sqrt(1.4*Pe_/rhoe_))^2)  % modified by gamma correction in Aftosmis Skylon Plumes

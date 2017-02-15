%Lift Force interpolator

%this module takes values from communicator and communicator-trim and finds
%appropriate AoA and flap deflection values to equalise a given q
%normalised lift force. this allows for splines to be created that only
%require M and lift force to give flap deflection, AoA and total drag
clear all;		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SPARTAN_SCALE
SPARTAN_SCALE = 1 % volumetric scale

global scattered
% Counts iterations of DIDO
global iterative_V
iterative_V = [];
global iterative_t
iterative_t = [];
global iteration
iteration = 1;
global iterative_V_f
iterative_V_f = [];
% Copy the current setting to archive
% This saves the entire problem file every time the program is run. 
% Disable this for conservatin of hard drive space
Timestamp = datestr(now,30)
copyfile('SecondStageProb.m',sprintf('../ArchivedResults/SecondStageProb_%s.m',Timestamp))
copyfile('SecondStageCost.m',sprintf('../ArchivedResults/SecondStageCost_%s.m',Timestamp))


% Inputs ============================================
%Take inputs of communicator matrices, these should be .txt files 
communicator = importdata('communicator.txt');
communicator_trim = importdata('communicator_trim.txt');


scattered.flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
scattered.flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
scattered.flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,6));


[MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';

scattered.Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
scattered.Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
scattered.pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');

% Produce Atmosphere Data
global Atmosphere
Atmosphere = dlmread('atmosphere.txt');

%produce scattered interpolants for thrust and fuel usage
% enginedata = dlmread('engineoutput_matrix');

%Produce scattered interpolants for vehicle data

%%
%NEWMODEL
shockdata = dlmread('ShockMat');
[MList,AOAList] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
pres_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
temp_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
scattered.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
scattered.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
scattered.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');
%%



scattered.data = dlmread('RESTM12DATA.txt');  
data = scattered.data;

newdata = [];
j=1;

% this sets the interpolation region! VERY IMPORTANT
% The interpolators have trouble with equivalence ratio because its 1 over
% a certain M
% this makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)
for i = 1: length(data(:,1))
    if data(i,1) < 5.
        newdata(j,:) = data(i,:);
        j=j+1;
    end
end

% IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,6));
% IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,3));

 p=polyfitn([data(:,1),data(:,2)],data(:,3),3)

M_englist = unique(sort(data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = floor(M_englist(1)):0.1:ceil(M_englist(end)); % enlarge spread, this is not necessary if you have a lot of engine data
M_eng_interp = unique(sort(data(:,1)));

T_englist = unique(sort(data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = floor(T_englist(1)):10:ceil(T_englist(end)); 
T_eng_interp = unique(sort(data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);
% grid.Isp_eng = IspScattered(grid.Mgrid_eng,grid.T_eng); % An 'interpolator' which only interpolates at the data points. This is just an easy way to make a grid.

for i = 1:30
    for j= 1:30
grid.Isp_eng(i,j) = polyvaln(p,[grid.Mgrid_eng(i,j) grid.T_eng(i,j)]);
    end
end

scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','spline');

scattered.equivalence = scatteredInterpolant(newdata(:,1),newdata(:,2),newdata(:,4), 'linear');
grid.eq_eng = scattered.equivalence(grid.Mgrid_eng,grid.T_eng);
scattered.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');




disp('Calculating Flight Dynamics Regime');

flapdeflection_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,3));
flapdrag_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,5));
flaplift_spline = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,4),communicator_trim(:,6));


[MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
Cl_Grid = reshape(communicator(:,3),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
Cd_Grid = reshape(communicator(:,4),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
pitchingmoment_Grid = reshape(communicator(:,11),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';

Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'spline','linear');
Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'spline','linear');
pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'spline','linear');


A = 62.77*SPARTAN_SCALE^(2/3); % reference area (m^2)

% searchs method, to search for a variety of M and lift forces
%equalises the pitching moment of flap and body to calculate lift. works
%towards the correct AoA (and corresponding flap pitching moment)
liftarray = [];
j = 1;
% for v = 1500:10:3000 % Velocity (m/s)
for v = 1500:25:3000 % Velocity (m/s)
todisp = [num2str(j/length(1500:25:3000)*100),' % complete '];
disp(todisp)
j = j+1;
%     for alt = 22000:100:40000 % Altitude (m)
for alt = 20000:250:50000 % Altitude (m)

%         for Lift = 0:2500:200000 % Lift force (N)   max mass of vehicle is 8755.1
Lift = 25000:2500:150000; % Lift force (N)   max mass of vehicle is 8755.1

Alphatemp = [];
flapdeflection=[];
Drag=[];
body_pitchingmoment=[];
error=[];
  
parfor i=1:length(Lift)

            
            c = spline( Atmosphere(:,1),  Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data

            rho = spline( Atmosphere(:,1),  Atmosphere(:,4), alt); % Calculate density using atmospheric data

            q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

            M = v./c; % Calculating Mach No (Descaled)
            
            %% Calculate Thrust Component ==================================
            
            %% Determine AoA ==============================================================
            T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), alt); 
            P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), alt); 

           
            Alphatemp(i) = fminsearch(@(Alpha)LiftError(M, Alpha, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift(i),T0,P0),5);

            %Fuel Cost ===========================================================================
            
            [Isp,Fueldt] = RESTM12int(M, Alphatemp(i), scattered, SPARTAN_SCALE,T0,P0);

            Thrust = Isp.*Fueldt*9.81;

            %======================================================================

            Cl1 = Cl_spline(M,Alphatemp(i));

            body_pitchingmoment(i) = pitchingmoment_spline(M, Alphatemp(i));% first approximation of pitchingmoment using only body lift

            Flap_lift = q./50000*flaplift_spline(M,Alphatemp(i),-body_pitchingmoment(i))*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments

            total_lift = Cl1*A*q + Flap_lift + Thrust*sin(deg2rad(Alphatemp(i))); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

            error(i) = abs(total_lift - Lift(i));
            

            flapdeflection(i) = flapdeflection_spline(M,Alphatemp(i),-body_pitchingmoment(i));

            Drag(i) = Cd_spline(M,Alphatemp(i))*A*q +  q/50000*flapdrag_spline(M,Alphatemp(i),-body_pitchingmoment(i))*SPARTAN_SCALE^(2/3);


end 

            liftarray = [liftarray;[v*ones(length(Lift),1),alt*ones(length(Lift),1),Lift.',Alphatemp.',flapdeflection.',Drag.',-body_pitchingmoment.',error.']];
    end
end

dlmwrite('liftarray',liftarray)


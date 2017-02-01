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

% =========================================================================
% SET RUN MODE
% =========================================================================
% Change const to set the target of the simulation. Much of the problem
% definition changes with const.

% const = 1x: No end constraint, used for optimal trajectory calculation
% const = 1: 50kPa limit, 12: 55 kPa limit, 13: 45 kPa limit, 14: 50kPa limit & 10% additional drag

% const = 3: Fuel mass is constrained at end point, used for constant
% dynamic pressure calculation (50kPa constrained)
% const = 31: simple model for guess calc 

global const
const = 1

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
%OLDMODEL
% [MList,AOAList] = ndgrid(unique(communicator(:,1)),unique(communicator(:,2)));
% M1_Grid = reshape(communicator(:,12),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% pres_Grid = reshape(communicator(:,13),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% temp_Grid = reshape(communicator(:,14),[length(unique(communicator(:,2))),length(unique(communicator(:,1)))]).';
% scattered.M1gridded = griddedInterpolant(MList,AOAList,M1_Grid,'spline','linear');
% scattered.presgridded = griddedInterpolant(MList,AOAList,pres_Grid,'spline','linear');
% scattered.tempgridded = griddedInterpolant(MList,AOAList,temp_Grid,'spline','linear');
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
% IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,6));
IspScattered = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
M_englist = unique(sort(data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = floor(M_englist(1)):0.1:ceil(M_englist(end)); % enlarge spread, this is not necessary if you have a lot of engine data
% M_eng_interp = unique(sort(data(:,1)));

T_englist = unique(sort(data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = floor(T_englist(1)):10:ceil(T_englist(end)); 
% T_eng_interp = unique(sort(data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);
grid.Isp_eng = IspScattered(grid.Mgrid_eng,grid.T_eng);
scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'spline','linear');
% scattered.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.Isp_eng,'linear','linear');
% scattered.equivalence = scatteredInterpolant(data(:,1),data(:,2),data(:,7));
scattered.equivalence = scatteredInterpolant(data(:,1),data(:,2),data(:,4));
% grid.eq_eng = scattered.equivalence(grid.Mgrid_eng,grid.T_eng);
% scattered.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'spline','linear');




disp('Calculating Flight Dynamics Regime');

%read communicator and create interpolatior / extrapolator splines
% Cl_spline = scatteredInterpolant(communicator(:,1),communicator(:,2),communicator(:,3)); %find cl given M, AoA
% % flaplift_spline_2 = scatteredInterpolant(communicator_trim(:,1),communicator_trim(:,2),communicator_trim(:,3),communicator_trim(:,6)); %find flap lift, given M, AoA, flap deflection
% Cd_spline = scatteredInterpolant(communicator(:,1),communicator(:,2),communicator(:,4)); % find Cd given M, AoA 
% % Alpha_spline = scatteredInterpolant(communicator(:,1),communicator(:,3),communicator(:,2)); % find AoA given M, Cl
% pitchingmoment_spline = scatteredInterpolant(communicator(:,1),communicator(:,2),communicator(:,11)); % find pitching moment given M, AoA changed from main code
% % pitchingmoment_spline = scatteredInterpolant(communicator(:,1),communicator(:,3),communicator(:,6));
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

% Cl_spline = griddedInterpolant(MList,AOAList,Cl_Grid,'linear','linear');
% Cd_spline = griddedInterpolant(MList,AOAList,Cd_Grid,'linear','linear');
% pitchingmoment_spline = griddedInterpolant(MList,AOAList,pitchingmoment_Grid,'linear','linear');


A = 62.77*SPARTAN_SCALE^(2/3); % reference area (m^2)

% golden sections method, to search for a variety of M and lift forces
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
Efficiency = [];
            
            c = spline( Atmosphere(:,1),  Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data

            rho = spline( Atmosphere(:,1),  Atmosphere(:,4), alt); % Calculate density using atmospheric data

            q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

            M = v./c; % Calculating Mach No (Descaled)
            
            %% Calculate Thrust Component ==================================
            
            if const == 1 || const == 14


                Efficiency = zeros(1,length(q));
                for j = 1:length(q)
                    if q(j) < 50000
                    Efficiency(j) = rho/(50000*2/v^2); % dont change this

                    else
                    Efficiency(j) = 1; % for 50kPa

                    end
                end
            elseif const == 12

                Efficiency = zeros(1,length(q));
                for j = 1:length(q)
                    if q(j) < 55000
                    Efficiency(j) = rho/(50000*2/v^2); % dont change this

                    else
                    Efficiency(j) = 1.1; % for 55kPa

                    end
                end
            elseif const == 13

                Efficiency = zeros(1,length(q));
                for j = 1:length(q)
                    if q(j) < 45000
                        Efficiency(j) = rho/(50000*2/v^2); % dont change this

                    else
                        Efficiency(j) = .9; % for 45kPa

                    end
                end
            elseif const == 3 || const == 31

            Efficiency = rho./(50000*2./v.^2); % linear rho efficiency, scaled to rho at 50000kpa

            end

Efficiency = 1;

            %% Determine AoA ==============================================================
            T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), alt); 
            P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), alt); 

           
            Alphatemp(i) = fminsearch(@(Alpha)LiftError(M, Alpha, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift(i),T0,P0),5);

%             Alpha = fminbnd(@(Alpha)LiftError(M, Alpha, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift,T0,P0),2,10);

            %             error = LiftError(M, Alpha, t_ratio, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline);
            

    %         Alpha = Alpha_spline(M, Liftq/A); % first approximation of alpha using only body lift

            %Fuel Cost ===========================================================================
            
            [Isp,Fueldt] = RESTM12int(M, Alphatemp(i), Efficiency, scattered, SPARTAN_SCALE,T0,P0);

            Thrust = Isp.*Fueldt*9.81;

            %======================================================================

            Cl1 = Cl_spline(M,Alphatemp(i));

            body_pitchingmoment(i) = pitchingmoment_spline(M, Alphatemp(i));% first approximation of pitchingmoment using only body lift

            Flap_lift = q./50000*flaplift_spline(M,Alphatemp(i),-body_pitchingmoment(i))*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments

            total_lift = Cl1*A*q + Flap_lift + Thrust*sin(deg2rad(Alphatemp(i))); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

            error(i) = abs(total_lift - Lift(i));
            
%             if error(i)>1000
%              v
%              alt
%              Lift
%              Alpha(i)
%             end
            

            flapdeflection(i) = flapdeflection_spline(M,Alphatemp(i),-body_pitchingmoment(i));

            Drag(i) = Cd_spline(M,Alphatemp(i))*A*q +  q/50000*flapdrag_spline(M,Alphatemp(i),-body_pitchingmoment(i))*SPARTAN_SCALE^(2/3);
           
            
%               Drag = Cd_spline(M,Alpha4)*A*q ; % changed to just body drag



end 
%             liftarray(end+1,1) = v;
%             liftarray(end,2) = alt;
%             liftarray(end,3) = Lift;
% liftarray(end,4) = Alpha;
% 
%             liftarray(end,5) = flapdeflection;
% 
%             liftarray(end,6) = Drag;
%             
%             liftarray(end,7) = -body_pitchingmoment;
%             
%             liftarray(end,8) = error;
            liftarray = [liftarray;[v*ones(length(Lift),1),alt*ones(length(Lift),1),Lift.',Alphatemp.',flapdeflection.',Drag.',-body_pitchingmoment.',error.']];
    end
end

dlmwrite('liftarray',liftarray)
% create splines
% given M and lift force / q , find AoA, flap deflection and total drag
% force / q
% if const == 3
% AoA_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,4)); 
% flapdeflection_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,5));
% Drag_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,6));
% Flap_pitchingmoment_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,7));
% % 
% else

% [vList,altList,liftList] = ndgrid(unique(liftarray(:,1)),unique(liftarray(:,2)),unique(liftarray(:,3)));
% 
% AoA_Grid = permute(reshape(liftarray(:,4),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
% AoA_spline = griddedInterpolant(vList,altList,liftList,AoA_Grid,'spline','linear');
% 
% flapdeflection_Grid = permute(reshape(liftarray(:,5),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
% flapdeflection_spline = griddedInterpolant(vList,altList,liftList,flapdeflection_Grid,'spline','linear');
% 
% Drag_Grid = permute(reshape(liftarray(:,6),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
% Drag_spline = griddedInterpolant(vList,altList,liftList,Drag_Grid,'spline','linear');
% 
% Flap_pitchingmoment_Grid = permute(reshape(liftarray(:,7),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
% Flap_pitchingmoment_spline = griddedInterpolant(vList,altList,liftList,Flap_pitchingmoment_Grid,'spline','linear');

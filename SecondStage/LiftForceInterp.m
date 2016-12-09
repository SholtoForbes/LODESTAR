function [AoA_spline, flapdeflection_spline, Drag_spline,Flap_pitchingmoment_spline] = LiftForceInterp(communicator,communicator_trim,const, Atmosphere, scattered, SPARTAN_SCALE)
%Lift Force interpolator

%this module takes values from communicator and communicator-trim and finds
%appropriate AoA and flap deflection values to equalise a given q
%normalised lift force. this allows for splines to be created that only
%require M and lift force to give flap deflection, AoA and total drag

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


A = 62.77; % reference area (m^2)

% golden sections method, to search for a variety of M and lift forces
%equalises the pitching moment of flap and body to calculate lift. works
%towards the correct AoA (and corresponding flap pitching moment)
liftarray = [];
j = 1;
for v = 1500:100:3000 % Velocity (m/s)
% for v = [1500 1700:25:3000] % Velocity (m/s)
    
todisp = [num2str(j/length(1500:100:3000)*100),' % complete '];
disp(todisp)
j = j+1;
    for alt = 20000:1000:50000 % Altitude (m)
%     for alt = 25000:250:40000 % Altitude (m)
        for Lift = 50000:5000:200000 % Lift force (N)   max mass of vehicle is 8755.1
%         for Lift = 0:2500:200000 % Lift force (N)   max mass of vehicle is 8755.1
            liftarray(end+1,1) = v;
            liftarray(end,2) = alt;
            liftarray(end,3) = Lift;
            
            c = spline( Atmosphere(:,1),  Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data

            rho = spline( Atmosphere(:,1),  Atmosphere(:,4), alt); % Calculate density using atmospheric data

            q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

            M = v./c; % Calculating Mach No (Descaled)
            
            %% Calculate Thrust Component ==================================
            
            kpa50_alt = interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000./v.^2);% altitude at each velocity, if it were a constant 50kPa q trajectory (used to compare with communicator matrix results) 
            kpa50_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), kpa50_alt); % Calculate density using atmospheric data
            
            if const == 1 || const == 14

                temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), alt);
                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 50000
                    Efficiency(i) = rho/(50000*2/v^2); % dont change this
                    t_ratio(i) = temp_actual(i)/kpa50_temp(i);
                    else
                    Efficiency(i) = 1; % for 50kPa
                    t_ratio(i) = 1; % 
                    end
                end
            elseif const == 12
                constq_alt = interp1(Atmosphere(:,4),Atmosphere(:,1),2*55000./v.^2); % altitude at each velocity, if it were a constant q trajectory (used to compare with communicator matrix results) 
                constq_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), constq_alt); % Calculate density using atmospheric data
                temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), alt);
                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 55000
                    Efficiency(i) = rho/(50000*2/v^2); % dont change this
                    t_ratio(i) = temp_actual(i)/kpa50_temp(i);
                    else
                    Efficiency(i) = 1.1; % for 55kPa
                    t_ratio(i) = constq_temp(i)/kpa50_temp(i);
                    end
                end
            elseif const == 13
                constq_alt = interp1(Atmosphere(:,4),Atmosphere(:,1),2*45000./v.^2); % altitude at each velocity, if it were a constant q trajectory (used to compare with communicator matrix results) 
                constq_temp =  spline( Atmosphere(:,1),  Atmosphere(:,2), constq_alt); % Calculate density using atmospheric data
                temp_actual = spline( Atmosphere(:,1),  Atmosphere(:,2), alt);
                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 45000
                        Efficiency(i) = rho/(50000*2/v^2); % dont change this
                        t_ratio(i) = temp_actual(i)/kpa50_temp(i);
                    else
                        Efficiency(i) = .9; % for 45kPa
                        t_ratio(i) = constq_temp(i)/kpa50_temp(i);
                    end
                end
            elseif const == 3 || const == 31
            Efficiency = rho./(50000*2./v.^2); % linear rho efficiency, scaled to rho at 50000kpa
            t_ratio = temp_actual/constq_temp;
            end



            %% Determine AoA ==============================================================
            
            
            

    %         Alpha = Alpha_spline(M, Liftq/A); % first approximation of alpha using only body lift

            Alpha1 = 0;
            
            %Fuel Cost ===========================================================================

%             Fueldt = FuelF_spline(M,Alpha1).*Efficiency;
% 
%             Isp = ThrustF_spline(M,Alpha1)./FuelF_spline(M,Alpha1); % this isnt quite Isp (doesnt have g) but doesnt matter
% 
%             Thrust = Isp.*Fueldt; % Thrust (N)
            

[Isp,Fueldt] = RESTM12int(M, Alpha1, t_ratio, Efficiency, scattered, SPARTAN_SCALE);

Thrust = Isp.*Fueldt*9.81;

% Thrust =  gridded.T_eng(M,Alpha1).*cos(deg2rad(Alpha1)).*Efficiency;
% Fueldt =  gridded.fuel_eng(M,Alpha1).*Efficiency;
            %======================================================================

            Cl1 = Cl_spline(M,Alpha1);

            body_pitchingmoment1 = pitchingmoment_spline(M, Alpha1);% first approximation of pitchingmoment using only body lift

            Flap_lift1 = q./50000*flaplift_spline(M,Alpha1,-body_pitchingmoment1);% first approximation of flap lift

            total_lift1 = Cl1*A*q + Flap_lift1 + Thrust*sin(deg2rad(Alpha1)); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq




            Alpha2 = 10; %first guesses of AoA
            %Fuel Cost ===========================================================================

%             Fueldt = FuelF_spline(M,Alpha2).*Efficiency;
% 
%             Isp = ThrustF_spline(M,Alpha2)./FuelF_spline(M,Alpha2); % this isnt quite Isp (doesnt have g) but doesnt matter
% 
%             Thrust = Isp.*Fueldt; % Thrust (N)



[Isp,Fueldt] = RESTM12int(M,Alpha2, t_ratio, Efficiency, scattered, SPARTAN_SCALE);

Thrust = Isp.*Fueldt*9.81;

% Thrust =  gridded.T_eng(M,Alpha2).*cos(deg2rad(Alpha2)).*Efficiency;
% Fueldt =  gridded.fuel_eng(M,Alpha2).*Efficiency;
            %======================================================================

            Cl2 = Cl_spline(M,Alpha2);

            body_pitchingmoment2 = pitchingmoment_spline(M, Alpha2);% first approximation of pitchingmoment using only body lift

            Flap_lift2 = q./50000*flaplift_spline(M,Alpha2,-body_pitchingmoment2);% first approximation of flap lift

            total_lift2 = Cl2*A*q + Flap_lift2 + Thrust*sin(deg2rad(Alpha2));



            Alpha3 = Alpha2 - (-1+sqrt(5))/2*(Alpha2-Alpha1); %first golden section point

            Alpha4 = Alpha1 + (-1+sqrt(5))/2*(Alpha2-Alpha1); %first guesses of AoA

            total_lift4 = 10000;
            total_lift3 = 0;

%             while abs(Alpha4 - Alpha3) > 0.05
%              while   abs(abs(Lift - total_lift4) - abs(Lift - total_lift3)) > 1
                 
             while   abs(Lift - total_lift4) > 100
            %Fuel Cost ===========================================================================

%             Fueldt = FuelF_spline(M,Alpha3).*Efficiency;
% 
%             Isp = ThrustF_spline(M,Alpha3)./FuelF_spline(M,Alpha3); % this isnt quite Isp (doesnt have g) but doesnt matter
% 
%             Thrust = Isp.*Fueldt; % Thrust (N)



[Isp,Fueldt] = RESTM12int(M,Alpha3, t_ratio, Efficiency, scattered, SPARTAN_SCALE);

Thrust = Isp.*Fueldt*9.81;

% Thrust =  gridded.T_eng(M,Alpha3).*cos(deg2rad(Alpha3)).*Efficiency;
% Fueldt =  gridded.fuel_eng(M,Alpha3).*Efficiency;
            %======================================================================
            Cl3 = Cl_spline(M,Alpha3);

            body_pitchingmoment3 = pitchingmoment_spline(M, Alpha3);% first approximation of pitchingmoment using only body lift

            Flap_lift3 = q./50000*flaplift_spline(M,Alpha3,-body_pitchingmoment3);% first approximation of flap lift

            total_lift3 = Cl3*A*q + Flap_lift3 + Thrust*sin(deg2rad(Alpha3));


            %Fuel Cost ===========================================================================

%             Fueldt = FuelF_spline(M,Alpha4).*Efficiency;
% 
%             Isp = ThrustF_spline(M,Alpha4)./FuelF_spline(M,Alpha4); % this isnt quite Isp (doesnt have g) but doesnt matter
% 
%             Thrust = Isp.*Fueldt; % Thrust (N)



[Isp,Fueldt] = RESTM12int(M,Alpha4, t_ratio, Efficiency, scattered, SPARTAN_SCALE);

Thrust = Isp.*Fueldt*9.81;

% Thrust =  gridded.T_eng(M,Alpha4).*cos(deg2rad(Alpha4)).*Efficiency;
% Fueldt =  gridded.fuel_eng(M,Alpha4).*Efficiency;
            %======================================================================
            Cl4 = Cl_spline(M,Alpha4);

            body_pitchingmoment4 = pitchingmoment_spline(M, Alpha4);% first approximation of pitchingmoment using only body lift

            Flap_lift4 = q./50000*flaplift_spline(M,Alpha4,-body_pitchingmoment4);% first approximation of flap lift

            total_lift4 = Cl4*A*q + Flap_lift4 + Thrust*sin(deg2rad(Alpha4));



                if abs(Lift - total_lift4) > abs(Lift - total_lift3)
                    Alpha2 = Alpha4;
                    Alpha4 = Alpha3;
                    Alpha3 = Alpha2 - (-1+sqrt(5))/2*(Alpha2-Alpha1);
                else
                    Alpha1 = Alpha3;
                    Alpha3 = Alpha4;
                    Alpha4 = Alpha1 + (-1+sqrt(5))/2*(Alpha2-Alpha1);
                end
                Alpha4
                M
                Alpha3

abs(Lift - total_lift4)
            end
            abs(Lift - total_lift4)

            flapdeflection = flapdeflection_spline(M,Alpha4,-body_pitchingmoment4);

            Drag = Cd_spline(M,Alpha4)*A*q +  q/50000*flapdrag_spline(M,Alpha4,-body_pitchingmoment4);
           
            
%               Drag = Cd_spline(M,Alpha4)*A*q ; % changed to just body drag

            liftarray(end,4) = Alpha4;

            liftarray(end,5) = flapdeflection;

            liftarray(end,6) = Drag;
            
            liftarray(end,7) = -body_pitchingmoment4;

        end
    end
end

% create splines
% given M and lift force / q , find AoA, flap deflection and total drag
% force / q

AoA_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,4)); 
flapdeflection_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,5));
Drag_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,6));
Flap_pitchingmoment_spline = scatteredInterpolant(liftarray(:,1),liftarray(:,2),liftarray(:,3),liftarray(:,7));


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
% end
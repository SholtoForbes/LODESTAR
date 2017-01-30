function [AoA_spline, flapdeflection_spline, Drag_spline,Flap_pitchingmoment_spline,liftarray] = LiftForceInterp(communicator,communicator_trim,const, Atmosphere, scattered, SPARTAN_SCALE)
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


A = 62.77*SPARTAN_SCALE^(2/3); % reference area (m^2)

% golden sections method, to search for a variety of M and lift forces
%equalises the pitching moment of flap and body to calculate lift. works
%towards the correct AoA (and corresponding flap pitching moment)
liftarray = [];
j = 1;
for v = 1500:100:3000 % Velocity (m/s)
% for v = 1500:50:3000 % Velocity (m/s)
todisp = [num2str(j/length(1500:100:3000)*100),' % complete '];
disp(todisp)
j = j+1;
    for alt = 20000:1000:50000 % Altitude (m)
% for alt = 20000:500:50000 % Altitude (m)

        for Lift = 0:5000:200000 % Lift force (N)   max mass of vehicle is 8755.1
% for Lift = 0:2500:200000 % Lift force (N)   max mass of vehicle is 8755.1

            liftarray(end+1,1) = v;
            liftarray(end,2) = alt;
            liftarray(end,3) = Lift;
            
            c = spline( Atmosphere(:,1),  Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data

            rho = spline( Atmosphere(:,1),  Atmosphere(:,4), alt); % Calculate density using atmospheric data

            q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

            M = v./c; % Calculating Mach No (Descaled)
            
            %% Calculate Thrust Component ==================================
            
            if const == 1 || const == 14


                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 50000
                    Efficiency(i) = rho/(50000*2/v^2); % dont change this

                    else
                    Efficiency(i) = 1; % for 50kPa

                    end
                end
            elseif const == 12

                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 55000
                    Efficiency(i) = rho/(50000*2/v^2); % dont change this

                    else
                    Efficiency(i) = 1.1; % for 55kPa

                    end
                end
            elseif const == 13

                Efficiency = zeros(1,length(q));
                for i = 1:length(q)
                    if q(i) < 45000
                        Efficiency(i) = rho/(50000*2/v^2); % dont change this

                    else
                        Efficiency(i) = .9; % for 45kPa

                    end
                end
            elseif const == 3 || const == 31

            Efficiency = rho./(50000*2./v.^2); % linear rho efficiency, scaled to rho at 50000kpa

            end



            %% Determine AoA ==============================================================
            T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), alt); 
            P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), alt); 
            
            Alpha = fminsearch(@(Alpha)LiftError(M, Alpha, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift,T0,P0),5);

%             Alpha = fminbnd(@(Alpha)LiftError(M, Alpha, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift,T0,P0),2,10);

            %             error = LiftError(M, Alpha, t_ratio, Efficiency, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline);
            

    %         Alpha = Alpha_spline(M, Liftq/A); % first approximation of alpha using only body lift

            %Fuel Cost ===========================================================================
            
            [Isp,Fueldt] = RESTM12int(M, Alpha, Efficiency, scattered, SPARTAN_SCALE,T0,P0);

            Thrust = Isp.*Fueldt*9.81;

            %======================================================================

            Cl1 = Cl_spline(M,Alpha);

            body_pitchingmoment = pitchingmoment_spline(M, Alpha);% first approximation of pitchingmoment using only body lift

            Flap_lift = q./50000*flaplift_spline(M,Alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments

            total_lift = Cl1*A*q + Flap_lift + Thrust*sin(deg2rad(Alpha)); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

            error = abs(total_lift - Lift);
            
            if error>1000
             v
             alt
             Lift
             Alpha
            end
            

            flapdeflection = flapdeflection_spline(M,Alpha,-body_pitchingmoment);

            Drag = Cd_spline(M,Alpha)*A*q +  q/50000*flapdrag_spline(M,Alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);
           
            
%               Drag = Cd_spline(M,Alpha4)*A*q ; % changed to just body drag

            liftarray(end,4) = Alpha;

            liftarray(end,5) = flapdeflection;

            liftarray(end,6) = Drag;
            
            liftarray(end,7) = -body_pitchingmoment;
            
            liftarray(end,8) = error;

        end
    end
end

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

[vList,altList,liftList] = ndgrid(unique(liftarray(:,1)),unique(liftarray(:,2)),unique(liftarray(:,3)));

AoA_Grid = permute(reshape(liftarray(:,4),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
AoA_spline = griddedInterpolant(vList,altList,liftList,AoA_Grid,'spline','linear');

flapdeflection_Grid = permute(reshape(liftarray(:,5),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
flapdeflection_spline = griddedInterpolant(vList,altList,liftList,flapdeflection_Grid,'spline','linear');

Drag_Grid = permute(reshape(liftarray(:,6),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
Drag_spline = griddedInterpolant(vList,altList,liftList,Drag_Grid,'spline','linear');

Flap_pitchingmoment_Grid = permute(reshape(liftarray(:,7),[length(unique(liftarray(:,3))),length(unique(liftarray(:,2))),length(unique(liftarray(:,1)))]),[3 2 1]);
Flap_pitchingmoment_spline = griddedInterpolant(vList,altList,liftList,Flap_pitchingmoment_Grid,'spline','linear');
end